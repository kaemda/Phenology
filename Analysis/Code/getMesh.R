library(sdmTMB)
library(sdmTMBextra)
library(fmesher)
library(sf)
library(INLA)

# MAKE MESH--------------------------------------------------------------------
getMesh <- function(data, data.sf) {
  
  northAmerica <- st_read("Data/Shapefiles/NorthAmerica/boundary_p_v2.shp") %>%
    subset(COUNTRY != "water/agua/d'eau" & COUNTRY != "FN") %>% 
    group_by(COUNTRY) %>% 
    summarize(geometry = st_union(geometry)) %>%
    st_transform(5070)
  
  # Calculate max edge. The spatial range is about 1/3 of the study area;
  # the max.edge is about 1/5 of that. 
  spatial.range = diff(range(data$Y))/3
  max.edge = spatial.range/5
  
  # Calculate boundary around points
  bnd <- INLA::inla.nonconvex.hull(cbind(data$X, data$Y), convex = -0.1)
  
  # Create INLA mesh with inside and outside edges
  mesh_inla <- fm_mesh_2d_inla(
    boundary = bnd,
    max.edge = c(1,2) * max.edge,
    offset = c(1,2) * max.edge,
    cutoff = max.edge/5
  )
  
  # Run sdmTMB mesh creation function with INLA mesh
  mesh <- make_mesh(data, c("X", "Y"), mesh = mesh_inla)
  plot(mesh)
  
  # Add North America as a barrier
  if (st_crs(data.sf) == st_crs(northAmerica)) {
    mesh_barrier <- sdmTMBextra::add_barrier_mesh(spde_obj = mesh, barrier_sf = northAmerica, range_fraction = 0.1, proj_scaling = 1000, plot = T)
    
    # Plot barrier to verify accuracy
    mesh_df_water <- mesh_barrier$mesh_sf[mesh_barrier$normal_triangles, ]
    mesh_df_land <- mesh_barrier$mesh_sf[mesh_barrier$barrier_triangles, ]
    
    pdf(file = paste0("Figures/",data$scientific_name[1],"/barrier_mesh.pdf"), width = 5, height = 5,)
    print(ggplot(northAmerica) +
            geom_sf() +
            geom_sf(data = mesh_df_water, size = 1, colour = "blue") +
            geom_sf(data = mesh_df_land, size = 1, colour = "green"))
    dev.off()
    
    return(mesh_barrier)
  } else {
    print("CRS of dataset and barrier are not the same!")
  }
  
  return(mesh_barrier)
  
}
