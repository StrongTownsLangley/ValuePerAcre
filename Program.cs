using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Newtonsoft.Json;
using Newtonsoft.Json.Linq;

/* STRONG TOWNS LANGLEY VALUE-PER-ACRE TOOL
 * Coded by: James Hansen (james@strongtownslangley.org)
 * https://github.com/StrongTownsLangley/ValuePerAcre
 * strongtownslangley.org
 */

namespace vpa
{
    class Program
    {
        // Arguments
        public static string mode = "";
        public static string taxRatesFile = null;
        public static string assessmentsFile = null;
        public static string assessmentsFilePIDField = "PID";
        public static string valuesFile = null;
        public static string parcelsFile = null;
        
        public static string outputFolder = "json";
        public static int taxRateDivider = 1000;
        public static int levels = 50;
        public static double blockSize = 100;
     
        // Specifications
        private static double minLongitude = 10000000;
        private static double maxLongitude = -10000000;
        private static double minLatitude = 10000000;
        private static double maxLatitude = -10000000;
        private static int yBlocks;
        private static int xBlocks;

        // Shared Classes
        class Property // For parcel based system
        {
            public string PID { get; set; } //do nulls as zero value
            public List<JToken> Features { get; set; } // change to list            
            public double Value { get; set; }
            public double TotalArea { get; set; }
            public double ValuePerArea { get; set; }
            public int Level { get; set; }
        };

        // Shared Functions
        private static double Distance(double lat1, double lon1, double lat2, double lon2, string unit)
        {
            double theta = lon1 - lon2;
            double dist = Math.Sin(Deg2Rad(lat1)) * Math.Sin(Deg2Rad(lat2)) + Math.Cos(Deg2Rad(lat1)) * Math.Cos(Deg2Rad(lat2)) * Math.Cos(Deg2Rad(theta));
            dist = Math.Acos(dist);
            dist = Rad2Deg(dist);
            double miles = dist * 60 * 1.1515;
            unit = unit.ToUpper();

            if (unit == "K")
                return miles * 1.609344;
            else if (unit == "N")
                return miles * 0.8684;
            else
                return miles;
        }
        public static double lastPercent = 0;
        private static void ProgressBar(int index, int totalIterations)
        {
            // Calculate the percentage progress
            double progress = ((double)index + 1) / totalIterations * 100;
            double percent = Math.Round(progress);

            if (lastPercent != percent)
            {
                // Update the progress bar
                int barLength = (int)Math.Round(progress / 2); // 50 characters for 100%
                Console.Write("[" + new string('=', barLength) + new string(' ', 50 - barLength) + "]");
                Console.Write(" " + percent + "%");
                Console.Write("\r"); // Move the cursor to the beginning of the line
                lastPercent = percent;
            }
        }

        // Constants for WGS84
        private const double a = 6378137.0; // Semi-major axis
        private const double e = 0.081819191; // First eccentricity
        public static double[] UtmToLatLon(double easting, double northing, int zone, bool isNorthernHemisphere)
        {
            if (!isNorthernHemisphere)
            {
                northing -= 10000000.0;
            }

            double k0 = 0.9996;
            double eccSquared = e * e;
            double eccPrimeSquared = eccSquared / (1 - eccSquared);
            double e1 = (1 - Math.Sqrt(1 - eccSquared)) / (1 + Math.Sqrt(1 - eccSquared));
            double N1, T1, C1, R1, D, M;
            double LongOrigin = (zone - 1) * 6 - 180 + 3;  // +3 puts origin in middle of zone
            double mu, phi1Rad;

            // Remove 500,000 meter offset for longitude
            double x = easting - 500000.0;
            double y = northing;

            // Calculate the Meridional Arc
            M = y / k0;

            // Calculate the footprint latitude
            mu = M / (a * (1 - eccSquared / 4 - 3 * eccSquared * eccSquared / 64 - 5 * eccSquared * eccSquared * eccSquared / 256));

            phi1Rad = mu + (3 * e1 / 2 - 27 * e1 * e1 * e1 / 32) * Math.Sin(2 * mu)
                         + (21 * e1 * e1 / 16 - 55 * e1 * e1 * e1 * e1 / 32) * Math.Sin(4 * mu)
                         + (151 * e1 * e1 * e1 / 96) * Math.Sin(6 * mu);

            double phi1 = phi1Rad * 180.0 / Math.PI;

            N1 = a / Math.Sqrt(1 - eccSquared * Math.Sin(phi1Rad) * Math.Sin(phi1Rad));
            T1 = Math.Tan(phi1Rad) * Math.Tan(phi1Rad);
            C1 = eccPrimeSquared * Math.Cos(phi1Rad) * Math.Cos(phi1Rad);
            R1 = a * (1 - eccSquared) / Math.Pow(1 - eccSquared * Math.Sin(phi1Rad) * Math.Sin(phi1Rad), 1.5);
            D = x / (N1 * k0);

            double lat = phi1Rad - (N1 * Math.Tan(phi1Rad) / R1) * (D * D / 2 - (5 + 3 * T1 + 10 * C1 - 4 * C1 * C1 - 9 * eccPrimeSquared) * D * D * D * D / 24
                    + (61 + 90 * T1 + 298 * C1 + 45 * T1 * T1 - 252 * eccPrimeSquared - 3 * C1 * C1) * D * D * D * D * D * D / 720);
            lat = lat * 180.0 / Math.PI;

            double lng = (D - (1 + 2 * T1 + C1) * D * D * D / 6
                 + (5 - 2 * C1 + 28 * T1 - 3 * C1 * C1 + 8 * eccPrimeSquared + 24 * T1 * T1) * D * D * D * D * D / 120) / Math.Cos(phi1Rad);
            lng = LongOrigin + lng * 180.0 / Math.PI;

            return new double[] { lng, lat };
        }
        // Function to populate blocks by level
        public static Dictionary<int, List<Block>> PopulateBlocksByLevel(Dictionary<int, Dictionary<int, double>> data, Dictionary<string, double> specArray, double valueCap, double degLat, double degLon, int numLevels)
        {
            var blocksByLevel = Enumerable.Range(0, numLevels).ToDictionary(i => i, _ => new List<Block>());

            foreach (var x in Enumerable.Range(0, (int)specArray["width_in_blocks"]))
            {
                foreach (var y in Enumerable.Range(0, (int)specArray["height_in_blocks"]))
                {
                    if (!data.ContainsKey(x) || !data[x].ContainsKey(y))
                        continue;
                    var value = data[x][y];
                    value = Math.Min(value, valueCap); // Limit value to valueCap

                    int level = 0;
                    if (value > 0)
                        level = Math.Min((int)Math.Floor(value / valueCap * numLevels), numLevels - 1); // Calculate level

                    // Convert block size from meters to degrees
                    double lonDeg = specArray["min_longitude"] + (x * (specArray["size_of_block"] / degLon));
                    double latDeg = specArray["max_latitude"] - (y * (specArray["size_of_block"] / degLat));

                    // Calculate latitude and longitude for each corner of the block
                    double minLon = lonDeg;
                    double maxLon = lonDeg + (specArray["size_of_block"] / degLon);
                    double minLat = latDeg;
                    double maxLat = latDeg + (specArray["size_of_block"] / degLat);

                    // Store block in corresponding level                    
                    blocksByLevel[level].Add(new Block( minLon, maxLon, minLat, maxLat, data[x][y] ) );
                }
            }

            return blocksByLevel;
        }

        private static double Deg2Rad(double deg)
        {
            return deg * Math.PI / 180.0;
        }

        private static double Rad2Deg(double rad)
        {
            return rad * 180.0 / Math.PI;
        }

        public class Block
        {
            public Block(double MinLon, double MaxLon, double MinLat, double MaxLat, double Value)
            {
                this.MinLon = MinLon;
                this.MaxLon = MaxLon;
                this.MinLat = MinLat;
                this.MaxLat = MaxLat;

                if (Values == null)
                    Values = new List<double>();
                Values.Add(Value);
            }

            public double MinLon { get; set; }
            public double MaxLon { get; set; }
            public double MinLat { get; set; }
            public double MaxLat { get; set; }
            public List<double> Values { get; set; }
        }


        // Main Program
        public static void Main(string[] args)
        {
            Console.WriteLine("Strong Towns Langley Value-per-Acre Tool");
            Console.WriteLine("https://github.com/StrongTownsLangley/ValuePerAcre");
            Console.WriteLine("https://strongtownslangley.org");
            Console.WriteLine("Coded by James Hansen (james@strongtownslangley.org)");
            Console.WriteLine("---------------------------------------------------");

            // Read Arguments //
            #region Read Arguments
            foreach (var arg in args)
            {
                if (arg.StartsWith("-from-tax-rates="))
                {
                    taxRatesFile = arg.Substring("-from-tax-rates=".Length);
                }
                else if (arg.StartsWith("-from-assessments="))
                {
                    assessmentsFile = arg.Substring("-from-assessments=".Length);
                }
                else if (arg.StartsWith("-from-values="))
                {
                    valuesFile = arg.Substring("-from-values=".Length);
                }
                else if (arg.StartsWith("-parcels="))
                {
                    parcelsFile = arg.Substring("-parcels=".Length);
                }
                else if (arg.StartsWith("-assessment-pid-field="))
                {
                    assessmentsFilePIDField = arg.Substring("-assessment-pid-field=".Length); // Default = "PID"
                }
                else if (arg.StartsWith("-output-folder="))
                {
                    outputFolder = arg.Substring("-output-folder=".Length);
                }
                else if (arg.StartsWith("-tax-rate-divider="))
                {
                    int.TryParse(arg.Substring("-tax-rate-divider=".Length), out taxRateDivider);
                }
                else if (arg.StartsWith("-levels="))
                {
                    int.TryParse(arg.Substring("-levels=".Length), out levels);
                }
                else if (arg.StartsWith("-block-size="))
                {
                    double.TryParse(arg.Substring("-block-size=".Length), out blockSize);
                }
            }

            if (!string.IsNullOrEmpty(taxRatesFile) && !string.IsNullOrEmpty(assessmentsFile))
            {
                mode = "assessments";

                if (!File.Exists(taxRatesFile))
                {
                    Console.WriteLine("Error: Tax Rates File \"" + taxRatesFile + "\" does not exist.");
                    return;
                }

                if (!File.Exists(assessmentsFile))
                {
                    Console.WriteLine("Error: Assessment File \"" + assessmentsFile + "\" does not exist.");
                    return;
                }
            }

            if (!string.IsNullOrEmpty(valuesFile))
            {
                mode = "values";

                if (!File.Exists(valuesFile))
                {
                    Console.WriteLine("Error: Values File \"" + valuesFile + "\" does not exist.");
                    return;
                }
            }

            if (mode == "")
            {
                Console.WriteLine("(1) Usage from Property Assessments (Have to Calculate Tax using Tax Rates File):");
                Console.WriteLine("    vpa.exe -from-tax-rates=\"tax-rates.json\" -from-assessments=\"assessment-file.geojson\" [-assessment-pid-field=\"PID\"] [-parcels=\"parcel-file.geojson\"]");
                Console.WriteLine("(2) Usage from Values (Tax or Value already calculated):");
                Console.WriteLine("    vpa.exe -from-values=\"value-file.geojson\" [-parcels=\"parcel-file.geojson\"]");                
                Console.WriteLine(@"-- Expected Values JSON Format (PID field only required if using parcel-file)
[
  {
    ""Latitude"": 49.00444377963,
    ""Longitude"": -122.64299023049,
    ""Value"": 577.68
    ""PID"": ""008-401-123""
  },
  {
    ""Latitude"": 49.00445057939,
    ""Longitude"": -122.63949034915,
    ""Value"": 8156.32
    ""PID"": ""008-401-456""
  },
  ...
]");
                Console.WriteLine("Optional Flags: [output-folder=\"json\"] [-tax-rate-divider=1000] [-levels=50] [-block-size=100]");
                return;
            }





            if (!Directory.Exists(outputFolder))
            {
                Directory.CreateDirectory(outputFolder);                
            }


            
            Console.WriteLine("Output Folder: {0}", outputFolder);            
            Console.WriteLine("Block Size: {0} m^2", blockSize);
            Console.WriteLine("Number of Levels to Group Blocks Into: {0}", levels);


            #endregion

            JArray valuesJArray = new JArray();

            // Calculate Tax //            
            #region Calculate Tax (if Usage 1) or Load Values (if Usage 2)
            if (mode == "assessments")
            {
                Console.WriteLine("From Tax File: {0}", taxRatesFile);
                Console.WriteLine("From Assessment File: {0}", assessmentsFile);
                Console.WriteLine("Tax Rate Divider: {0}", taxRateDivider);

                string taxRateStr = File.ReadAllText(taxRatesFile);
                Dictionary<string, double> taxRate = JsonConvert.DeserializeObject<Dictionary<string, double>>(taxRateStr);

                // Load JSON
                JObject assessmentsJObject = null;
                using (FileStream fs = new FileStream(assessmentsFile, FileMode.Open))
                using (StreamReader sr = new StreamReader(fs))
                using (JsonTextReader reader = new JsonTextReader(sr))
                {
                    assessmentsJObject = JObject.Load(reader);
                }

                // Output File


                var list = assessmentsJObject["features"];
                foreach (var entry in list)
                {
                    var properties = entry["properties"];

                    // Calculate tax
                    double tax = 0;
                    foreach (var taxType in taxRate.Keys)
                    {
                        tax += Convert.ToDouble(properties[taxType + "_Buildings"]) * (taxRate[taxType] / taxRateDivider);
                        tax += Convert.ToDouble(properties[taxType + "_Land"]) * (taxRate[taxType] / taxRateDivider);
                    }
                    tax = Math.Round(tax, 2);

                    //Console.WriteLine($"PID {properties["PID"]} - Address: {properties["Unit"]} {properties["House"]} {properties["Street"]} - Lot Size: {properties["Lot_Size"]} - Tax amount: ${tax}");

                    dynamic dataEntry = new JObject();
                    //dataEntry.PID = properties["PID"];
                    //dataEntry.Lot_Size = properties["Lot_Size"];
                    //dataEntry.Lot_Desc = properties["Lot_Desc"];

                    dataEntry.Latitude = properties["Latitude"];
                    dataEntry.Longitude = properties["Longitude"];
                    dataEntry.Value = tax;
                    if(properties["PID"] != null)
                        dataEntry.PID = properties["PID"];

                    valuesJArray.Add(dataEntry);
                }

                valuesFile = Path.GetFileNameWithoutExtension(assessmentsFile) + ".Values.geojson";
                Console.WriteLine("Writing Values JSON '" + valuesFile + "'...");
                File.WriteAllText(Path.Combine(outputFolder, valuesFile), JsonConvert.SerializeObject(valuesJArray, Formatting.Indented));
            }
            else
            {
                Console.WriteLine("From Values File: {0}", valuesFile);

                // Load JSON                
                using (FileStream fs = new FileStream(valuesFile, FileMode.Open))
                using (StreamReader sr = new StreamReader(fs))
                using (JsonTextReader reader = new JsonTextReader(sr))
                {
                    valuesJArray = JArray.Load(reader);
                }

            }
            #endregion

            var geoJsonArray = new List<object>();

            // Data Output
            Dictionary<string, double> groupedByBlockSpecArray = new Dictionary<string, double>();
            foreach (JObject entry in valuesJArray)
            {
                double longitude = (double)entry["Longitude"];
                double latitude = (double)entry["Latitude"];

                if (longitude == 0 || latitude == 0)
                    continue;

                minLongitude = Math.Min(minLongitude, longitude);
                maxLongitude = Math.Max(maxLongitude, longitude);
                minLatitude = Math.Min(minLatitude, latitude);
                maxLatitude = Math.Max(maxLatitude, latitude);
            }

            Console.WriteLine("Min Longitude: " + minLongitude);
            Console.WriteLine("Max Longitude: " + maxLongitude);
            Console.WriteLine("Min Latitude: " + minLatitude);
            Console.WriteLine("Max Latitude: " + maxLatitude);

            groupedByBlockSpecArray["min_longitude"] = minLongitude;
            groupedByBlockSpecArray["max_longitude"] = maxLongitude;
            groupedByBlockSpecArray["min_latitude"] = minLatitude;
            groupedByBlockSpecArray["max_latitude"] = maxLatitude;

            double mWidth = Distance(minLatitude, minLongitude, maxLatitude, minLongitude, "K") * 1000;
            double mHeight = Distance(minLatitude, minLongitude, minLatitude, maxLongitude, "K") * 1000;

            Console.WriteLine("Width: " + mWidth + " metres");
            Console.WriteLine("Height: " + mHeight + " metres");

            groupedByBlockSpecArray["width_in_metres"] = mWidth;
            groupedByBlockSpecArray["height_in_metres"] = mHeight;

            if (parcelsFile == null)
            {
                Console.WriteLine("Method: Points-based mode (No Parcel File Specified)");
                #region Calculate Into Blocks
             


                yBlocks = (int)Math.Round(mWidth / blockSize);
                xBlocks = (int)Math.Round(mHeight / blockSize);

                Console.WriteLine("Width: " + xBlocks + " " + blockSize + "m^2 blocks");
                Console.WriteLine("Height: " + yBlocks + " " + blockSize + "m^2 blocks");

                groupedByBlockSpecArray["size_of_block"] = blockSize;
                groupedByBlockSpecArray["width_in_blocks"] = xBlocks;
                groupedByBlockSpecArray["height_in_blocks"] = yBlocks;

                // Calculate tax values for blocks
                Dictionary<int, Dictionary<int, double>> groupedByBlockArray = new Dictionary<int, Dictionary<int, double>>();
                foreach (var entry in valuesJArray)
                {
                    double longitude = Convert.ToDouble(entry["Longitude"]);
                    double latitude = Convert.ToDouble(entry["Latitude"]);
                    double tax = Convert.ToDouble(entry["Value"]);

                    int yBlockPos = yBlocks - (int)Math.Round(Distance(minLatitude, minLongitude, latitude, minLongitude, "K") * 1000 / blockSize);
                    int xBlockPos = (int)Math.Round(Distance(minLatitude, minLongitude, minLatitude, longitude, "K") * 1000 / blockSize);

                    if (!groupedByBlockArray.ContainsKey(xBlockPos))
                        groupedByBlockArray[xBlockPos] = new Dictionary<int, double>();
                    if (!groupedByBlockArray[xBlockPos].ContainsKey(yBlockPos))
                        groupedByBlockArray[xBlockPos][yBlockPos] = 0;
                    groupedByBlockArray[xBlockPos][yBlockPos] += tax;
                }

                string groupedByBlockFile = Path.GetFileNameWithoutExtension(valuesFile) + ".ValuesByBlockData.geojson";
                string groupedByBlockSpecFile = Path.GetFileNameWithoutExtension(valuesFile) + ".ValuesByBlockSpec.geojson";

                Console.WriteLine("Writing Value By Block Data JSON '" + groupedByBlockFile + "'...");
                Console.WriteLine("Writing Value By Block Spec JSON '" + groupedByBlockSpecFile + "'...");

                File.WriteAllText(Path.Combine(outputFolder, groupedByBlockFile), JsonConvert.SerializeObject(groupedByBlockArray, Formatting.Indented));
                File.WriteAllText(Path.Combine(outputFolder, groupedByBlockSpecFile), JsonConvert.SerializeObject(groupedByBlockSpecArray, Formatting.Indented));
                #endregion

                #region Build GeoJSON
                // Build GeoJSON Files grouped by Level //
                // Calculate size of one degree of latitude and longitude in meters
                double degLat = 111320; // meters
                double degLon = 111320 * Math.Cos(Math.PI * (groupedByBlockSpecArray["min_latitude"] + groupedByBlockSpecArray["max_latitude"]) / 360); // meters

                // Loop through all the blocks sorted from highest value to lowest value
                var sortedValues = groupedByBlockArray.Values.SelectMany(dict => dict.Values).OrderByDescending(value => value).ToList();

                // Variables to store the best distribution and corresponding block value
                Dictionary<int, int> bestDistribution = null;
                double bestValue = 0;
                int totalIterations = sortedValues.Count();
                int index = 0;
                double lastPercent = 0;
                // Initialize the progress bar
                Console.WriteLine("Searching for best distribution...");
                foreach (var testBlockValue in sortedValues)
                {
                    // Calculate the percentage progress
                    ProgressBar(index, totalIterations);


                    // Populate blocks by level for the current test block value
                    var blocksByLevel = PopulateBlocksByLevel(groupedByBlockArray, groupedByBlockSpecArray, testBlockValue, degLat, degLon, levels);

                    // Calculate how evenly distributed the values are in each level
                    var distribution = blocksByLevel.ToDictionary(kv => kv.Key, kv => kv.Value.Count);

                    // Check if this distribution is better than the previous best distribution
                    if (bestDistribution == null || distribution.Values.Max() - distribution.Values.Min() < bestDistribution.Values.Max() - bestDistribution.Values.Min())
                    {
                        bestDistribution = distribution;
                        bestValue = testBlockValue;
                    }

                    index++;
                }
                Console.WriteLine(""); // Clear Progress Bar

                // Use the best block value to create the final GeoJSON output
                var finalBlocksByLevel = PopulateBlocksByLevel(groupedByBlockArray, groupedByBlockSpecArray, bestValue, degLat, degLon, levels);
                geoJsonArray = new List<object>();
                var levelsInfoArray = new List<object>();
                // Generate GeoJSON for each level
                foreach (var level in finalBlocksByLevel.Keys)
                {
                    var minLevelBlockValue = finalBlocksByLevel[level].Min(m => m.Values.Sum());
                    var maxLevelBlockValue = finalBlocksByLevel[level].Max(m => m.Values.Sum());
                    var avgLevelBlockValue = finalBlocksByLevel[level].Average(m => m.Values.Sum());
                    var numLevelBlockValue = finalBlocksByLevel[level].Sum(m => m.Values.Count);
                    var geoJson = new
                    {
                        type = "FeatureCollection",
                        features = new List<object>(),
                        info = new List<object>()
                    };

                    // Loop through all blocks in this level
                    foreach (var block in finalBlocksByLevel[level])
                    {
                        var minBlockValue = block.Values.Min();
                        var maxBlockValue = block.Values.Max();
                        var avgBlockValue = block.Values.Average();
                        var numBlockValue = block.Values.Count;

                        var feature = new
                        {
                            type = "Feature",
                            properties = new { level = level, minBlockValue = minBlockValue, maxBlockValue = maxBlockValue, avgBlockValue = avgBlockValue, numBlockValue = numBlockValue },
                            geometry = new
                            {
                                type = "Polygon",
                                coordinates = new[]
                                {
                                    new[]
                                    {
                                        new[] { block.MinLon, block.MinLat },
                                        new[] { block.MaxLon, block.MinLat },
                                        new[] { block.MaxLon, block.MaxLat },
                                        new[] { block.MinLon, block.MaxLat },
                                        new[] { block.MinLon, block.MinLat }
                                    }
                                }
                            }
                        };
                        geoJson.features.Add(feature);
                    }

                    // Store info about this level
                    var info = new { level = level, minLevelBlockValue = minLevelBlockValue, maxLevelBlockValue = maxLevelBlockValue, avgLevelBlockValue = avgLevelBlockValue, numLevelBlockValue = numLevelBlockValue };
                    geoJson.info.Add(info);
                    levelsInfoArray.Add(info);

                    // Store in Array
                    geoJsonArray.Add(geoJson);

                    // Write GeoJSON to file
                    var levelFile = "level_" + level + ".json";
                    Console.WriteLine("Writing Level JSON '" + levelFile + "' file...");
                    File.WriteAllText(Path.Combine(outputFolder, levelFile), Newtonsoft.Json.JsonConvert.SerializeObject(geoJson));

                }
                var levelInfoFile = "level_info.json";
                Console.WriteLine("Writing Level JSON '" + levelInfoFile + "' file...");
                File.WriteAllText(Path.Combine(outputFolder, levelInfoFile), Newtonsoft.Json.JsonConvert.SerializeObject(levelsInfoArray));
                #endregion
            }
            else
            {
                Console.WriteLine("Method: Parcel-based mode (Parcel File Specified)");
                #region Load Parcel File
                // Load the GeoJSON file
                Console.WriteLine("Reading Parcels File: {0}", parcelsFile);
                JObject parcelsJson = JObject.Parse(File.ReadAllText(parcelsFile));

                JArray features = (JArray)parcelsJson["features"];

                

                // Populate Properties Array with All Polygons for Property and their Total Area and Fix Polygons for Leaflet
                Console.WriteLine("Processing Polygons...");
                int index = 0;

                var PIDToValuesDictionary = valuesJArray
                .Where(f => f["PID"] != null)
                .GroupBy(f => f["PID"].ToString())
                .ToDictionary(g => g.Key, g => (double)g.First()["Value"]);
                Dictionary<string, Property> PIDToPropertiesDictionary = new Dictionary<string, Property>();


                foreach (var feature in features)
                {
                    ProgressBar(index, features.Count());
                    index++;

                    string PID = feature["properties"]["PID"].ToString();
                    var geometry = feature["geometry"];
                    double area = 0; 
                    double value = 0;
                    if (geometry["type"].ToString() == "Polygon")
                    {
                        var coordinates = geometry["coordinates"];
                        foreach (var polygon in coordinates)
                        {
                                                       
                            for (int i = 0; i < polygon.Count(); i++)
                            {
                                // Fix Coordinates into standard Lat/Long Format    
                                double[] thisCoord = ((JArray)polygon[i]).Select(x => (double)x).ToArray();
                                double[] transformedCoord = UtmToLatLon(thisCoord[0], thisCoord[1], 10, true);  // Assuming UTM zone 10N
                                polygon[i] = new JArray(transformedCoord[0], transformedCoord[1]);
                            }

                            for (int i = 0; i < polygon.Count(); i++)
                            {
                                // Calculate Area    
                                double[] thisCoord = ((JArray)polygon[i]).Select(x => (double)x).ToArray();
                                double[] nextCoord = ((JArray)polygon[(i + 1) % polygon.Count()]).Select(x => (double)x).ToArray();
                                                            
                                area += thisCoord[0] * nextCoord[1] - nextCoord[0] * thisCoord[1];
                            }
                            var areaInSquareDegrees = Math.Abs(area) / 2.0;
                            var areaInSquareMetres = areaInSquareDegrees * 12362500000;
                            area = areaInSquareMetres / 4046.85642; // Acres
                        }
                    }


                    if (PIDToValuesDictionary.ContainsKey(PID))
                        value = PIDToValuesDictionary[PID];

                    if (value == 0)
                        continue;

                    
                    if (!PIDToPropertiesDictionary.ContainsKey(PID))
                    {
                        var property = new Property()
                        {
                            PID = PID,
                            Features = new List<JToken>() { feature },
                            Value = value,
                            TotalArea = area
                        };
                        PIDToPropertiesDictionary.Add(PID,property);
                    } else {
                        PIDToPropertiesDictionary[PID].Features.Add(feature);
                        PIDToPropertiesDictionary[PID].TotalArea += area;
                    }                    
                }
                Console.WriteLine("");


                // Calculate Value Per Area
                foreach (var propertyKVP in PIDToPropertiesDictionary)
                {
                    var property = propertyKVP.Value;
                    property.ValuePerArea = property.Value / property.TotalArea;
                }

                // PID Key no longer needed, Convert to List
                var properties = PIDToPropertiesDictionary.Select(m => m.Value).OrderBy(x => x.ValuePerArea).ToList();

                //
                string valuePerAreaFile = Path.GetFileNameWithoutExtension(assessmentsFile) + ".ValuePerArea.geojson";
                Console.WriteLine("Writing ValuePerArea JSON '" + valuePerAreaFile + "'...");
                using (StreamWriter file = File.CreateText(Path.Combine(outputFolder, valuePerAreaFile)))
                using (JsonTextWriter writer = new JsonTextWriter(file))
                {
                    writer.Formatting = Formatting.Indented;
                    JsonSerializer serializer = new JsonSerializer();

                    writer.WriteStartArray();
                    foreach (var property in properties)
                    {
                        var propertyStripped = new Property();
                        propertyStripped.PID = property.PID;
                        propertyStripped.Value = property.Value;
                        propertyStripped.TotalArea = property.TotalArea;
                        propertyStripped.ValuePerArea = property.ValuePerArea;
                        propertyStripped.Level = property.Level;

                        serializer.Serialize(writer, propertyStripped);
                    }
                    writer.WriteEndArray();
                }

                Console.WriteLine("Finding best distribution...");

                    int totalItems = properties.Count;                    
                    int itemsPerLevel = totalItems / levels;
                    int remainder = totalItems % levels;

                    int currentIndex = 0;

                    for (int level = 0; level < levels; level++)
                    {
                        int itemsInThisLevel = itemsPerLevel + (level < remainder ? 1 : 0);

                        for (int i = 0; i < itemsInThisLevel; i++)
                        {
                            properties[currentIndex].Level = level;
                            currentIndex++;
                        }
                    }

                    // Create and populate GeoJSON objects for each group 
                    Dictionary<int, List<Property>> groupedPIDs = properties
                    .GroupBy(p => p.Level)
                    .ToDictionary(g => g.Key, g => g.ToList());
                    var levelsInfoArray = new List<object>();
                    foreach (var group in groupedPIDs)
                    {
                        int level = group.Key;
                        var levelProperties = group.Value;

                        var minLevelBlockValue = levelProperties.Min(m => m.ValuePerArea);
                        var maxLevelBlockValue = group.Value.Max(m => m.ValuePerArea);
                        var avgLevelBlockValue = group.Value.Average(m => m.ValuePerArea);
                        var numLevelBlockValue = group.Value.Count();
                        var geoJson = new
                        {
                            type = "FeatureCollection",
                            features = new List<object>(),
                            info = new List<object>()
                        };

                        foreach (var property in levelProperties)
                        {
                            if (property == null) continue;
                            foreach (var feature in property.Features)
                            {
                                feature["properties"] = new JObject
                                {
                                    { "level", level },
                                    { "value_per_acre", property.ValuePerArea }
                                };
                                geoJson.features.Add(feature);                                
                            }
                        }

                        // Store info about this level
                        var info = new { level = level, minLevelBlockValue = minLevelBlockValue, maxLevelBlockValue = maxLevelBlockValue, avgLevelBlockValue = avgLevelBlockValue, numLevelBlockValue = numLevelBlockValue };
                        geoJson.info.Add(info);
                        levelsInfoArray.Add(info);

                        // Write GeoJSON to file
                        var levelFile = "level_" + level + ".json";
                        Console.WriteLine("Writing Level JSON '" + levelFile + "' file...");
                        File.WriteAllText(Path.Combine(outputFolder, levelFile), Newtonsoft.Json.JsonConvert.SerializeObject(geoJson));
                        geoJsonArray.Add(geoJson);
                    }
                    var levelInfoFile = "level_info.json";
                    Console.WriteLine("Writing Level JSON '" + levelInfoFile + "' file...");
                    File.WriteAllText(Path.Combine(outputFolder, levelInfoFile), Newtonsoft.Json.JsonConvert.SerializeObject(levelsInfoArray));

                #endregion
            }
            Console.WriteLine();


            #region Write Map HTML
            // Replace {AVGLAT}, {AVGLON}, {LEVELS}, {DATALIST}
            // taxLevel.addData(data);
            var website = File.ReadAllText("vpa.template.html");

            website = website.Replace("{AVGLAT}", ((minLatitude + maxLatitude) / 2).ToString());
            website = website.Replace("{AVGLON}", ((minLongitude + maxLongitude) / 2).ToString());
            website = website.Replace("{LEVELS}", levels.ToString());
            

            string geoJsonString = "";
            var geoIndex = 0;
            foreach (var geoJson in geoJsonArray)
            {
                geoJsonString += "taxLevels[" + geoIndex + "].addData(" + Newtonsoft.Json.JsonConvert.SerializeObject(geoJson) + "); \r\n";
                geoIndex++;
            }

            string websiteStatic = website.Replace("{DATALIST}", geoJsonString);
            string websiteDynamic = website.Replace("{DATALIST}",
@"for (var i = 0; i < numLevels; i++) {
             fetch('level_' + i + '.json')
                .then(response => response.json())
                .then(data => taxLevel.addData(data))
                .catch(error => console.error('Error loading GeoJSON file:', error));
           
}"
);

            var websiteStaticFile = "website.static.html";
            Console.WriteLine("Writing Static Website Heatmap '" + websiteStaticFile + "' file...");
            File.WriteAllText(Path.Combine(outputFolder, websiteStaticFile), websiteStatic);

            var websiteDynamicFile = "website.dynamic.html";
            Console.WriteLine("Writing Dynamic Website Heatmap '" + websiteDynamicFile + "' file...");
            File.WriteAllText(Path.Combine(outputFolder, websiteDynamicFile), websiteDynamic);     
            #endregion

            Console.WriteLine("Complete!");
        }
    }
}
