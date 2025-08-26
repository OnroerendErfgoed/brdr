import osmnx as ox

# Stel het type infrastructuur in (bijv. 'cycleway' voor fietspaden)
tags = {"highway": "cycleway"}

# Download alle fietspaden in Leuven
leuven_fietspaden = ox.features_from_place("Leuven, Belgium", tags)

# Toon de eerste rijen
print(leuven_fietspaden.head())

# Optioneel: opslaan als GeoPackage of GeoJSON
leuven_fietspaden.to_file("leuven_fietspaden.gpkg", driver="GPKG")
