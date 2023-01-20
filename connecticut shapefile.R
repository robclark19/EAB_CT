require(rgdal)
require(ggplot2)

shp <- st_read("./Shape/Connecticut_Mainland_Polygon.shp", stringsAsFactors = F)

plot(shp)


shp <- st_read("./Shape/CT_County_Index.shp", stringsAsFactors = F)

st_geometry_type(shp)

plot(shp)
