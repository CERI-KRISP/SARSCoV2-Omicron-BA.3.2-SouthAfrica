# ---- packages ----
library(treeio)
library(ggtree)
library(dplyr)
library(ggplot2)
library(readr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ape)      # for branching.times if needed

## PATHS
#MCC_PATH   <- "../data/ba3_all/ba3.2/20251210/beast/ba32.country.gtr.ucln.skygrid.50.MCC.tre"
MCC_PATH   <- "../data/ba3_all/ba3.2/20251210/beast/travel informed/ba32.country.gtr.ucln.skygrid.50.MCC.tre" # travel informed
COORDS_CSV <- "../data/ba3_all/ba3.2/20251210/metadata/country_coords_corrected.csv" # travel informed
MRS_DECIMAL <- 2025.90958904 # Optional


## helpers functions
pick_first <- function(nms, choices) {
  cand <- intersect(choices, nms)
  if (length(cand) == 0) NA_character_ else cand[1]
}
normkey <- function(x) gsub("\\s+", " ", trimws(tolower(x)))

infer_mrsd_from_labels <- function(lbl) {
  yrs <- suppressWarnings(as.numeric(gsub("^.*?/(\\d{4})$", "\\1", lbl)))
  yrs <- yrs[!is.na(yrs)]
  if (length(yrs) == 0) NA_real_ else max(yrs) + 0.50
}

decimal_year_to_date <- function(x) {
  year  <- floor(x)
  frac  <- x - year
  as.Date(paste0(year, "-01-01")) + round(frac * 365.25)
}

## preprocessing
mcc <- read.beast(MCC_PATH)
df  <- as_tibble(mcc)

## choose string location column (TreeAnnotator put labels under `location.rate`)
pick_location <- function(x) {
  if ("location.rate" %in% names(x) && is.character(x$`location.rate`)) {
    x$`location.rate`
  } else if ("country" %in% names(x) && is.character(x$country)) {
    x$country
  } else if ("location" %in% names(x) && is.character(x$location)) {
    x$location
  } else {
    lbl <- if ("label" %in% names(x)) x$label else rep(NA_character_, nrow(x))
    sub(".*?/([^/]+)/[^/]+$", "\\1", lbl)
  }
}

df <- df %>%
  mutate(
    country_use = pick_location(cur_data_all()),
    is_tip      = !is.na(label)
  )

## get node heights
# Prefer BEAST-provided height_median; else height; else compute from phylo.
if ("height_median" %in% names(df)) {
  height_col <- "height_median"
} else if ("height" %in% names(df)) {
  height_col <- "height"
} else {
  height_col <- NULL
}

if (is.null(height_col) || all(is.na(df[[height_col]]))) {
  tr <- treeio::get.tree(mcc) # ape::phylo
  Ntip  <- length(tr$tip.label)
  Nnode <- tr$Nnode
  bt <- ape::branching.times(tr) # internal nodes ages before present
  node_heights <- numeric(Ntip + Nnode)
  node_heights[seq_len(Ntip)] <- 0
  node_heights[(Ntip+1):(Ntip+Nnode)] <- bt[as.character((Ntip+1):(Ntip+Nnode))]
  df$height_use <- node_heights[df$node]
} else {
  df$height_use <- as.numeric(df[[height_col]])
}

## derive MRSD
if (is.na(MRS_DECIMAL)) {
  tip_lbls <- df$label[df$is_tip]
  MRS_DECIMAL <- infer_mrsd_from_labels(tip_lbls)
  message(sprintf("MRSD inferred as %.2f from tip labels (edit MRS_DECIMAL to override).", MRS_DECIMAL))
}

## parent/child with transition times
df_edges <- df %>%
  filter(!is.na(parent)) %>%
  transmute(
    child_node     = node,
    parent_node    = parent,
    parent_country = country_use[match(parent, node)],
    child_country  = country_use,
    child_height   = height_use
  ) %>%
  filter(!is.na(parent_country), !is.na(child_country),
         parent_country != child_country) %>%
  mutate(
    child_time_year = MRS_DECIMAL - as.numeric(child_height)  # calendar year of change (approx.)
  )

## Tree
p_tree <- ggtree(mcc) %<+% df +
  geom_tree(aes(color = country_use), linewidth = 0.7) +
  theme_tree2(legend.position = "left") +
  labs(color = "Country")
print(p_tree)

#ggsave("../outputs/figure3/mcc_discrete_tree_20251210.pdf", height = 4, width = 4)
ggsave("../outputs/figure3/mcc_discrete_tree_travel_informed_20251210.pdf", height = 4, width = 4) # travel informed

## Phylogeograhpy
if (nrow(df_edges) > 0) {
  coords_raw <- read_csv(COORDS_CSV, show_col_types = FALSE)
  
  # detect coordinate columns automatically
  nms <- names(coords_raw)
  country_col <- pick_first(nms, c("country","location","name","Country","Location","Name"))
  lat_col     <- pick_first(nms, c("lat","latitude","Lat","Latitude","y","Y"))
  lon_col     <- pick_first(nms, c("lon","long","longitude","Lon","Long","Longitude","x","X"))
  stopifnot(!is.na(country_col), !is.na(lat_col), !is.na(lon_col))
  
  coords <- coords_raw %>%
    transmute(
      country = .data[[country_col]],
      lat     = as.numeric(.data[[lat_col]]),
      lon     = as.numeric(.data[[lon_col]])
    ) %>%
    filter(!is.na(country), !is.na(lat), !is.na(lon)) %>%
    distinct(country, .keep_all = TRUE) %>%
    mutate(key = normkey(country))
  
  moves <- df_edges %>%
    mutate(from_key = normkey(parent_country),
           to_key   = normkey(child_country)) %>%
    left_join(coords %>% select(from_key = key, lat_from = lat, lon_from = lon), by = "from_key") %>%
    left_join(coords %>% select(to_key   = key, lat_to   = lat, lon_to   = lon), by = "to_key") %>%
    filter(!is.na(lat_from), !is.na(lat_to), !is.na(lon_from), !is.na(lon_to))
  
  if (nrow(moves) == 0) {
    message("No movements could be joined to coordinates. Check country names in the CSV.")
  } else {
    world <- rnaturalearth::ne_countries(scale = "small", returnclass = "sf")
    
    # dynamic timescale
    rng   <- range(moves$child_time_year, na.rm = TRUE)
    pad   <- 0.05
    lims  <- c(rng[1] - pad, rng[2] + pad)
    # always label as Month Year
    labf <- function(x) format(decimal_year_to_date(x), "%b %Y")
    span  <- diff(lims)
    if (span <= 1) {
      brks <- seq(ceiling(lims[1]*10)/10, floor(lims[2]*10)/10, by = 0.1)
      labf <- scales::number_format(accuracy = 0.1)
    } else {
      brks <- seq(floor(lims[1]), ceiling(lims[2]), by = 1)
    }
    
    # colour palette
    pal <- rev(RColorBrewer::brewer.pal(9, "YlOrRd"))
    #pal <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
    
    ggplot() +
      geom_sf(data = world, fill = "grey95", color = "grey80", linewidth = 0.1) +
      geom_curve(
        data = moves,
        aes(x = lon_from, y = lat_from, xend = lon_to, yend = lat_to,
            color = child_time_year),
        curvature = 0.25, alpha = 0.85,
        linewidth = 0.9,
        arrow = arrow(type = "closed", length = unit(5, "pt"))
      ) +
      geom_point(data = coords, aes(x = lon, y = lat), size = 1) +
      scale_color_gradientn(
        name   = "Transition year",
        colors = pal,
        limits = lims,
        breaks = brks,
        labels = labf
      ) +
      coord_sf(crs = "EPSG:4326") +
      theme_minimal(base_size = 11) +
      theme(panel.grid = element_blank(), axis.title = element_blank())
  }
}

## save plot
#ggsave(filename = "../outputs/figure3/discrete_spread_map_20251210.png",
#  plot     = last_plot(), width    = 9, height = 5, dpi = 300)

#ggsave(filename = "../outputs/figure3/discrete_spread_map_20251210.pdf",
#  plot     = last_plot(), width    = 9, height = 5)

## save plot (travel informed)
ggsave(
  filename = "../outputs/figure3/discrete_spread_map_travel_informed_20251210.png",
  plot     = last_plot(),        # or 'p_map' if you assigned it
  width    = 9, height = 5, dpi = 300
)

ggsave(
  filename = "../outputs/figure3/discrete_spread_map_travel_informed_20251210.pdf",
  plot     = last_plot(),        # or 'p_map' if you assigned it
  width    = 9, height = 5
)


## Legend item


library(ggplot2)
library(RColorBrewer)
library(cowplot)  # for get_legend()

# build a dummy color scale using the same palette and breaks as your map
rng   <- range(moves$child_time_year, na.rm = TRUE)
pad   <- 0.05
lims  <- c(rng[1] - pad, rng[2] + pad)
span  <- diff(lims)

if (span <= 2) {
  brks <- seq(ceiling(lims[1]*2)/2, floor(lims[2]*2)/2, by = 0.5)
  labf <- function(x) format(decimal_year_to_date(x), "%b %Y")
} else {
  brks <- seq(floor(lims[1]), ceiling(lims[2]), by = 1)
  labf <- waiver()
}

pal <- rev(RColorBrewer::brewer.pal(9, "YlOrRd"))

# make a dummy plot that only has the color scale
legend_plot <- ggplot(data.frame(x = lims, y = lims, z = lims)) +
  geom_raster(aes(x, y, fill = z)) +
  scale_fill_gradientn(
    colours = pal,
    limits  = lims,
    breaks  = brks,
    labels  = labf,
    name    = "Transition date",
    guide   = guide_colorbar(
      direction = "horizontal",
      barwidth  = unit(5, "in"),
      barheight = unit(0.2, "in"),
      title.position = "top",
      title.hjust = 0.5
    )
  ) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  )

legend_plot

# extract and save only the legend
legend_only <- cowplot::get_legend(legend_plot)

ggsave(
  filename = "../outputs/figure3/discrete_spread_legend_20251210.pdf",
  plot     = legend_only,
  width    = 6, height = 1
)
