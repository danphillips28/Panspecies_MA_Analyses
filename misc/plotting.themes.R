# =========================
# Global plotting constants
# =========================
BASE_SIZE   <- 9
BASE_FAMILY <- "Arial"

LINE_SIZE   <- 0.6
POINT_SIZE  <- 1.8

# Okabe-Ito colourblind-safe palette
PALETTE_DISCRETE <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

# =========================
# Light theme
# =========================
theme_my_light <- function() {
  theme_minimal(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
    theme(
      # Axes
      axis.line = element_line(color = "black", linewidth = LINE_SIZE),
      axis.ticks = element_line(color = "black", linewidth = LINE_SIZE),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      scale_y_continuous(n.breaks = 10),
      
      # Panel
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = "black"),
      panel.grid.minor.y = element_line(color = "black", linewidth = 0.25),
      
      # Plot
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(size = BASE_SIZE, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = BASE_SIZE -1, face = "italic", hjust = 0),
      plot.caption = element_text(size = BASE_SIZE -1, hjust = 0),
      
      # Legend
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(color = "black"),
      legend.title = element_text(color = "black"),
      
      # Facets
      strip.background = element_rect(fill = "grey90", color = NA),
      strip.text = element_text(color = "black", face = "bold")
    )
}

# =========================
# Dark theme
# =========================
theme_my_dark <- function() {
  theme_minimal(base_size = BASE_SIZE, base_family = BASE_FAMILY) +
    theme(
      # Axes
      axis.line = element_line(color = "white", linewidth = LINE_SIZE),
      axis.ticks = element_line(color = "white", linewidth = LINE_SIZE),
      axis.text = element_text(color = "white"),
      axis.title = element_text(color = "white"),
      scale_y_continuous(n.breaks = 10),
     
      # Panel
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "black", color = NA),
      
      # Plot
      plot.background = element_rect(fill = "white", color = NA),
      plot.title = element_text(size = BASE_SIZE, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = BASE_SIZE -1, face = "italic", hjust = 0),
      plot.caption = element_text(size = BASE_SIZE -1, hjust = 0),
      
      # Legend
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.text = element_text(color = "white"),
      legend.title = element_text(color = "white"),
      
      # Facets
      strip.background = element_rect(fill = "grey20", color = NA),
      strip.text = element_text(color = "white", face = "bold")
    )
}

# =========================
# Scales
# =========================
scale_color_my <- function(...) {
  scale_color_manual(values = PALETTE_DISCRETE, ...)
}

scale_fill_my <- function(...) {
  scale_fill_manual(values = PALETTE_DISCRETE, ...)
}

# =========================
# Theme switcher
# =========================
set_my_theme <- function(mode = c("light", "dark")) {
  mode <- match.arg(mode)
  
  if (mode == "light") {
    theme_set(theme_my_light())
  } else {
    theme_set(theme_my_dark())
  }
  
  # Set consistent geom defaults
  update_geom_defaults("point", list(size = POINT_SIZE))
  update_geom_defaults("line", list(linewidth = LINE_SIZE))
}
