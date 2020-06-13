# PREP -------------------------------------------------------------------------
pacman::p_load(mapview, RColorBrewer, zoo, covid19.analytics, lwgeom, animation, magick,
               COVID19, data.table, PBSmapping, rnaturalearth, tidyverse,
               hablar, sf, feather, janitor, lubridate)
set_wd_to_script_path()
options(scipen = 9999)


# World map --------------------------------------------------------------------
# Read the data
# mini_world <- read_sf('data/ne_110m_land/ne_110m_land.shp')
mini_world <- ne_countries(scale = 50, returnclass = "sf")


# COVID data -------------------------------------------------------------------
d <- covid19.analytics::covid19.data("ts-deaths") %>%
  pivot_longer(c(-(1:4)), names_to = "date") %>%
  as_tibble() %>%
  convert(dte(date),
          chr(Country.Region)) %>%
  group_by(country = Country.Region, date) %>%
  summarise(deaths = sum_(value)) %>%
  group_by(country) %>%
  arrange(date) %>%
  mutate(deaths = deaths - lag(deaths),
         deaths_r7d = rollmean(deaths, 7, fill = NA, align = "right", na.rm = T)) %>%
  ungroup() %>%
  arrange(country, date) %>%
  drop_na(deaths_r7d) %>%
  mutate(country = case_when(
    country == "US"  ~ "United States of America",
    T ~ country
  ))

centers <- map_dfr(unique(d$date), function(.x) {
  sf <- mini_world %>%
    st_centroid() %>%
    select(name) %>%
    left_join(d %>% filter(date == .x), by = c("name" = "country")) %>%
    mutate(deaths_r7d = if_na(deaths_r7d, 0))

  sf <- sf %>%
    bind_cols(do.call(rbind, st_geometry(sf)) %>%
                as_tibble() %>% setNames(c("lon","lat"))) %>%
    summarise(center_lng = weighted.mean(lon, w = deaths_r7d),
              center_lat = weighted.mean(lat, w = deaths_r7d)) %>%
    st_drop_geometry() %>%
    mutate(date = .x)
})

# Plot function ----------------------------------------------------------------
# Much of the map projection in the following fuunctoin is copied from:
# https://gist.github.com/fzenoni/ef23faf6d1ada5e4a91c9ef23b0ba2c1

# Define the orthographic projection
# Choose lat_0 with -90 <= lat_0 <= 90 and lon_0 with -180 <= lon_0 <= 180

create_map <- function(.x) {
  print(paste("Total frames", length(frames)))
  print(paste("This frame no:", .x))
  lat <- smooth_df %>% filter(row == .x) %>% pull(center_lat)
  lon <- smooth_df %>% filter(row == .x) %>% pull(center_lng)
  .date <- smooth_df %>% filter(row == .x) %>% pull(date) %>% first()

  ortho <- paste0('+proj=ortho +lat_0=', lat, ' +lon_0=', lon, ' +x_0=0 +y_0=0 +a=6371000 +b=6371000 +units=m +no_defs')

  # Define the polygon that will help you finding the "blade"
  # to split what lies within and without your projection
  circle <- st_point(x = c(0,0)) %>% st_buffer(dist = 6371000) %>% st_sfc(crs = ortho)

  # Project this polygon in lat-lon
  circle_longlat <- circle %>% st_transform(crs = 4326)

  # circle_longlat cannot be used as it is
  # You must decompose it into a string with ordered longitudes
  # Then complete the polygon definition to cover the hemisphere
  if(lat != 0) {
    circle_longlat <- st_boundary(circle_longlat)

    circle_coords <- st_coordinates(circle_longlat)[, c(1,2)]
    circle_coords <- circle_coords[order(circle_coords[, 1]),]
    circle_coords <- circle_coords[!duplicated(circle_coords),]

    # Rebuild line
    circle_longlat <- st_linestring(circle_coords) %>% st_sfc(crs = 4326)

    if(lat > 0) {
      rectangle <- list(rbind(circle_coords,
                              c(X = 180, circle_coords[nrow(circle_coords), 'Y']),
                              c(X = 180, Y = 90),
                              c(X = -180, Y = 90),
                              c(X = -180, circle_coords[1, 'Y']),
                              circle_coords[1, c('X','Y')])) %>%
        st_polygon() %>% st_sfc(crs = 4326)
    } else {
      rectangle <- list(rbind(circle_coords,
                              c(X = 180, circle_coords[nrow(circle_coords), 'Y']),
                              c(X = 180, Y = -90),
                              c(X = -180, Y = -90),
                              c(X = -180, circle_coords[1, 'Y']),
                              circle_coords[1, c('X','Y')])) %>%
        st_polygon() %>% st_sfc(crs = 4326)
    }

    circle_longlat <- st_union(st_make_valid(circle_longlat), st_make_valid(rectangle))
  }

  # This visualization shows the visible emisphere in red
  # ggplot() +
  #   geom_sf(data = mini_world) +
  #   geom_sf(data = circle_longlat, color = 'red', fill = 'red', alpha = 0.3)

  # A small negative buffer is necessary to avoid polygons still disappearing in a few pathological cases
  # I should not change the shapes too much
  visible <- st_intersection(st_make_valid(mini_world), st_buffer(circle_longlat, -0.09)) %>%
    st_transform(crs = ortho)

  # DISCLAIMER: This section is the outcome of trial-and-error and I don't claim it is the best approach
  # Resulting polygons are often broken and they need to be fixed
  # Get reason why they're broken
  broken_reason <- st_is_valid(visible, reason = TRUE)

  # First fix NA's by decomposing them
  # Remove them from visible for now
  na_visible <- visible[is.na(broken_reason),]
  visible <- visible[!is.na(broken_reason),]

  # Open and close polygons
  na_visible <- st_cast(na_visible, 'MULTILINESTRING') %>%
    st_cast('LINESTRING', do_split=TRUE)
  na_visible <- na_visible %>% mutate(npts = npts(geometry, by_feature = TRUE))
  # Exclude polygons with less than 4 points
  na_visible <- na_visible %>%
    filter(npts >=4) %>%
    select(-npts) %>%
    st_cast('POLYGON')

  # Fix other broken polygons
  broken <- which(!st_is_valid(visible))
  for(land in broken) {
    result = tryCatch({
      # visible[land,] <- st_buffer(visible[land,], 0) # Sometimes useful sometimes not
      visible[land,] <- st_make_valid(visible[land,]) %>%
        st_collection_extract()
    }, error = function(e) {
      visible[land,] <<- st_buffer(visible[land,], 0)
    })
  }

  # Bind together the two tables
  visible <- rbind(visible, na_visible)


  # Filter out path of eye
  center_path <- smooth_df %>%
    st_as_sf(coords = c("center_lng", "center_lat"),
             crs = 4326) %>%
    filter(row <= .x)

  # FILL COVID
  fill_sf <- st_collection_extract(visible) %>%
    left_join(ddd %>% filter(row == .x), by = c("admin" = "country")) %>%
    mutate(deaths_r7d = if_na(deaths_r7d, 0),
           wght = deaths_r7d / pop_est * 1e6)

  # MAX fill
  max_fill <- st_collection_extract(visible) %>%
    left_join(ddd, by = c("admin" = "country")) %>%
    mutate(deaths_r7d = if_na(deaths_r7d, 0),
           wght = deaths_r7d / pop_est * 1e6) %>%
    pull(deaths_r7d) %>%
    max_()

  # Final plot
  ggplot() +
    geom_sf(data = circle,
            fill = 'steelblue', color = "black", size = .2) + 
    theme_void() +
    geom_sf(data = fill_sf,
            fill = "gray70", size = .12, color = "transparent") +
    geom_sf(data = fill_sf,
            aes(fill = deaths_r7d), color = "black", size = .12) +
    scale_fill_gradientn(colors = c("gray98", brewer.pal(9, "Reds")[1:9]), na.value = "gray98",
                         limits = c(0, max_fill), trans = "sqrt", guide = guide_colourbar(barwidth = unit(.4, "npc"),
                                                                                          barheight = unit(.03, "npc"),
                                                                                          title.position="top", title.hjust = 0.5)) +
    scale_alpha_continuous(range = c(.6, 1)) +
    geom_sf(data = center_path %>% st_coordinates() %>% st_linestring() %>% st_sfc(crs = 4326), color = "black", size = 1) +
    geom_sf(data = sf::st_point(x = c(lat, lon)) %>% st_sfc(crs = ortho), color = "black", fill = "#EF3B2C", size = 18, shape = 21, stroke = 1) +
    coord_sf(crs = ortho) +
    labs(title = "COVID-19: Center of gravity",
         subtitle = as.Date(.date) %>% format("%Y %B %d"),
         fill = "Daily deaths per country",
         caption = "\nSource: CSSE at Johns Hopkins University | Graphics: David Sjoberg\nNote 1: Center show the spatial weighted centroid for number of deaths (rolling 7 day avarage).\nNote 2: Red color indicate daily deaths per per country (square root rolling 7 day avarage).") +
    theme(legend.text = element_text(size = 25),
          legend.title = element_text(size = 25),
          plot.title = element_text(size = 90, hjust = .5, face = "bold"),
          plot.caption = element_text(size = 20, hjust = .5, color = "gray70"),
          plot.subtitle = element_text(size = 70, hjust = .5),
          legend.position = "bottom",
          plot.margin = margin(2, 2, 2, 2, "cm"))
}



# Smooth frames ----------------------------------------------------------------
.n_frames_smooth <- 8
smooth_df <- bind_rows(centers %>% slice(1, 1, 1),
                       map_dfr(2:nrow(centers), function(.r) {
                         .df <- centers %>%
                           slice((.r-1):.r)
                         tibble(center_lng = .df$center_lng[2] - (.df$center_lng[2] - .df$center_lng[1]) * seq(0, 1, length.out = .n_frames_smooth) %>% rev() %>% .[-1],
                                center_lat = .df$center_lat[2] - (.df$center_lat[2] - .df$center_lat[1]) * seq(0, 1, length.out = .n_frames_smooth) %>% rev() %>% .[-1],
                                date = .df$date[2])})) %>%
  mutate(row = row_number()) %>%
  as_tibble()

ddd <- map_dfr(d %>% split(d$country), function(.d) {
  bind_rows(.d %>% slice(1, 1, 1),
            map_dfr(2:nrow(.d), function(.r) {
              .df <- .d %>%
                slice((.r-1):.r)
              tibble(deaths = .df$deaths[2] - (.df$deaths[2] - .df$deaths[1]) * seq(0, 1, length.out = .n_frames_smooth) %>% rev() %>% .[-1],
                     deaths_r7d = .df$deaths_r7d[2] - (.df$deaths_r7d[2] - .df$deaths_r7d[1]) * seq(0, 1, length.out = .n_frames_smooth) %>% rev() %>% .[-1],
                     date = .df$date[2],
                     country = .df$country[1])})) %>%
    mutate(row = row_number()) %>%
    as_tibble()
})

# ** Test plot -----------------------------------------------------------------
# create_map(100)
# create_map(310)
# animation::saveVideo({
#   walk(c(340), ~create_map(.x) %>% print())
# }, movie.name = "mygif2.mp4", interval = 0.1, ani.width = 2000, ani.height = 2000)



# VIDEO ------------------------------------------------------------------------
# Be aware that it takes some time to render all plots. Maybe an hour.
frames <- c(1:nrow(smooth_df), rep(nrow(smooth_df), as.integer(nrow(smooth_df)/4)))
animation::saveVideo({
  walk(frames, ~create_map(.x) %>% print())
}, movie.name = "mygif2.mp4", interval = 0.04, ani.width = 2000, ani.height = 2000)
