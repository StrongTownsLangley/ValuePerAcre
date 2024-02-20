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
        public static string valuesFile = null;
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
                Console.WriteLine("(1) Usage from Property Assessments (Have to Calculate Tax using Tax Rates File): vpa.exe -from-tax-rates=\"tax-rates.json\" -from-assessments=\"assessment-file.geojson\"");
                Console.WriteLine("(2) Usage from Values (Tax or Value already calculated): vpa.exe -from-values=\"value-file.geojson\"");
                Console.WriteLine(@"-- Expected Values JSON Format:
[
  {
    ""Latitude"": 49.00444377963,
    ""Longitude"": -122.64299023049,
    ""Value"": 577.68
  },
  {
    ""Latitude"": 49.00445057939,
    ""Longitude"": -122.63949034915,
    ""Value"": 8156.32
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

            #region Calculate Into Blocks
            // Data Output
            Dictionary<string, double> groupedByBlockSpecArray = new Dictionary<string, double>();

            // Output File
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

            // Build GeoJSON Files grouped by Level //
            #region Build GeoJSON

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
            Console.Write("[                                                  ]\r"); // 50 characters for the progress bar

            foreach (var testBlockValue in sortedValues)
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
            Console.WriteLine("");

            // Use the best block value to create the final GeoJSON output
            var finalBlocksByLevel = PopulateBlocksByLevel(groupedByBlockArray, groupedByBlockSpecArray, bestValue, degLat, degLon, levels);
            var geoJsonArray = new List<object>();
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

            #region Write Map HTML
            // Replace {AVGLAT}, {AVGLON}, {LEVELS}, {DATALIST}
            // taxLevel.addData(data);
            var website = File.ReadAllText("vpa.template.html");

            website = website.Replace("{AVGLAT}", ((minLatitude + maxLatitude) / 2).ToString());
            website = website.Replace("{AVGLON}", ((minLongitude + maxLongitude) / 2).ToString());
            website = website.Replace("{LEVELS}", levels.ToString());
            

            string geoJsonString = "";
            index = 0;
            foreach (var geoJson in geoJsonArray)
            {
                geoJsonString += "taxLevels[" + index + "].addData(" + Newtonsoft.Json.JsonConvert.SerializeObject(geoJson) + "); \r\n";
                index++;
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
