# Value-per-Acre Mapping Tool
This tool takes a list of latitude and longitude point values, combined with either property assessments and tax rates or pre-calculated values, to produce an interactive web-based value-per-acre style heat-map where points and values are totalled into blocks.
![valueperacre_map](https://github.com/StrongTownsLangley/ValuePerAcre/assets/160652425/e6086d8a-2f75-4d5f-9fb8-9a3ffa4089b9)

## Demo
This tool was used to compile the map at https://strongtownslangley.org/maps?revenue-map

## Download
You can download the tool from the releases page: https://github.com/StrongTownsLangley/ValuePerAcre/releases/

## How to Use
This tool was designed to work with Township of Langley data specifically, but should be flexible enough to work for other places.

The tool has two modes, one where it will calculate the values from input tax rates and assessed property values, and another where it will take pre-calculated values.

```
(1) Usage from Property Assessments (Have to Calculate Tax using Tax Rates File): vpa.exe -from-tax-rates="tax-rates.json" -from-assessments="assessment-file.geojson"
(2) Usage from Values (Tax or Value already calculated): vpa.exe -from-values="value-file.geojson"
Optional Flags: [output-folder="json"] [-tax-rate-divider=1000] [-levels=50] [-block-size=100]"
```

In the **Output folder** it will create **level_\*.json** files each containing the data to be used as a layer in [Leaflet](https://github.com/Leaflet/Leaflet), as well viewable maps **website.static.html** (which contains all the layers data for local testing) and **website.dynamic.html** (which loads the level_\*.json files dynamically - NOTE: this will not work locally due to CORS) which load Leaflet and display the map.

### What is Value-per-Acre and Why?
Looking at a place where properties are grouped into blocks based on their taxable value is allows us to see what areas are the most and least productive and contribute the most and least revenue to city finances.

Value per Acre analysis is advocated by [Strong Towns](https://strongtowns.org) to evaluate the efficiency of land use by focusing on its productivity relative to infrastructure investment. By prioritizing compact, mixed-use development and identifying opportunities for infill and adaptive reuse, this can improve economic vitality, reduce sprawl, and create more sustainable, resilient communities. It aligns with Strong Towns' principles by encouraging long-term planning that improves the financial health of cities.

- Langley Township offers comprehensive assessed property values on its open data portal, which include latitude and longitude coordinates at *https://data-tol.opendata.arcgis.com/search?tags=Land%20and%20Parcel%20Information*
- The tool calculates taxable values by multiplying them with corresponding tax rates which are released in PDF format every year e.g. *https://www.tol.ca/en/services/resources/property-taxes/document-feed/2021-Tax-Rates.pdf*
- The total area is partitioned into 100m² blocks (which is adjustable using the -block-size flag)
- Each property value is added to its respective block based on it's coordinates.
- The tool searches for which property value to use as the highest value cap, the blocks are then grouped into 50 levels from highest to lowest (number of levels adjustable with -levels flag) according to the most even distribution.

**NOTE: This is a "point-based" mapping tool. As such parcel size and shape currently does not factor into the calculation. A large single site which pays a high amount of property tax may show as a single lone high-value block, which doesn't accurately show it's true value-per-acre spread out across multiple blocks. In future this will be improved to encorporate parcel sizes and shapes.**



### Tax Rates and Assessments JSON Format
If using method (1), then the tax rates and assessment file must be in the following format:

tax-rates.json:
```
{
    "Residential": 3.80248,
    "Utilities": 42.53820,
    "SupportiveHousing": 2.37008,
    "MajorIndustry": 8.34595,
    "LightIndustry": 10.24375,
    "Business": 12.02044,
    "ManagedForest": 0,
    "Rec_NonProfit": 6.73828,
    "Farm": 15.97448
}
```
These rates are then used in combination with the *_Land and *_Buildings values in the assessment file. Below is a version with just the required/used structure and fields:
assessment-file.geojson:
```
{
"type": "FeatureCollection",
"name": "Assessments",
"crs": { "type": "name", "properties": { "name": "" } },
"features": [
{ "type": "Feature", "properties": { "Residential_Buildings": 113700, "Residential_Land": 0, "Utilities_Improvements": 0, "Utilities_Land": 0, "SupportiveHousing_Buildings": 0, "SupportiveHousing_Land": 0, "MajorIndustry_Buildings": 0, "MajorIndustry_Land": 0, "LightIndustry_Buildings": 0, "LightIndustry_Land": 0, "Business_Buildings": 0, "Business_Land": 0, "ManagedForest_Improvements": 0, "ManagedForest_Land": 0, "Rec_NonProfit_Buildings": 0, "Rec_NonProfit_Land": 0, "Farm_Buildings": 0, "Farm_Land": 9098, "Latitude": 49.004443779630002, "Longitude": -122.64299023049 } },
{ "type": "Feature", "properties": { "Residential_Buildings": 191000, "Residential_Land": 1954000, "Utilities_Improvements": 0, "Utilities_Land": 0, "SupportiveHousing_Buildings": 0, "SupportiveHousing_Land": 0, "MajorIndustry_Buildings": 0, "MajorIndustry_Land": 0, "LightIndustry_Buildings": 0, "LightIndustry_Land": 0, "Business_Buildings": 0, "Business_Land": 0, "ManagedForest_Improvements": 0, "ManagedForest_Land": 0, "Rec_NonProfit_Buildings": 0, "Rec_NonProfit_Land": 0, "Farm_Buildings": 0, "Farm_Land": 0, "Latitude": 49.004450579390003, "Longitude": -122.63949034914999 } },
...
]
}
```

A simplifified version might look like this:
tax-rates.json:
```
{
    "Residential": 5
}
```
assessment-file.geojson:
```
{
"type": "FeatureCollection",
"name": "Assessments",
"crs": { "type": "name", "properties": { "name": "" } },
"features": [
{ "type": "Feature", "properties": { "Residential_Buildings": 113700, "Residential_Land": 0, "Latitude": 49.004443779630002, "Longitude": -122.64299023049 } },
{ "type": "Feature", "properties": { "Residential_Buildings": 191000, "Residential_Land": 1954000, "Latitude": 49.004450579390003, "Longitude": -122.63949034914999 } },
...
]
}
```

### Values JSON Format
If using method (2), the values file should be in this format. This is likely the method you will want to use if your data is from another city or source.
value-file.geojson:
```
[
  {
    "Latitude": 49.00444377963,
    "Longitude": -122.64299023049,
    "Value": 577.68
  },
  {
    "Latitude": 49.00445057939,
    "Longitude": -122.63949034915,
    "Value": 8156.32
  },
  ...
]
```

## Other Similar Projects
Fellow Strong Towns Local Conversation group [A Better Cobb](https://abettercobb.substack.com/) created their own set of Value per Acre tools, available at https://github.com/ABetterCobb/ValuePerAcre/

## Contributing

The project is a Visual Studio 2010 project in C# .NET 4.0. It's an old platform, but I like to use older platforms for simple tools like this.

If you encounter any issues while using it or have suggestions for improvement, please [open an issue](https://github.com/StrongTownsLangley/ValuePerAcre/issues) on the GitHub repository. Pull requests are also welcome.

## License

This program is released under the [Apache 2.0 License](https://github.com/StrongTownsLangley/ValuePerAcre/blob/main/LICENSE). If you use it for your website or project, please provide credit to **Strong Towns Langley** and preferably link to this GitHub.
