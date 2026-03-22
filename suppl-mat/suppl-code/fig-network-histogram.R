# plot interaction matrix and histograms
# generates a combined figure: networks (top) and coefficient densities (bottom)

# setup ------------------------------------------------------------------

# load setup
source("code/01-setup.R")

# load empirical matrices
load("data/processed/empirical/empirical-matrices.RData")


# positive and negative interaction networks -----------------------------

# create a combined adjacency matrix for the full network
A <- cbind(alpha, matrix(0, nrow(alpha), 6))
colnames(A) <- c(colnames(alpha), colnames(gamma))
B <- cbind(t(gamma), matrix(0, 6, 6))
colnames(B) <- c(colnames(t(gamma)), colnames(gamma))
AB <- as.matrix(rbind(A, B))

# assign node names for the full network (plants + herbivores)
node_names <- c(
  "Am",
  "Ao",
  "Ae",
  "Be",
  "Cj",
  "Ca",
  "Cr",
  "Dg",
  "Dc",
  "Er",
  "Em",
  "Fa",
  "Fr",
  "Gv",
  "Gd",
  "Gr",
  "Lv",
  "Lp",
  "Lc",
  "Ma",
  "Or",
  "Pe",
  "Ph",
  "Pl",
  "Pa",
  "Pp",
  "Pt",
  "Ra",
  "Rx",
  "Sp",
  "So",
  "To",
  "Tf",
  "Tp",
  "Vo",
  "Vp",
  "Gb",
  "Cd",
  "Ci",
  "Ee",
  "Pg",
  "Pp"
)
colnames(AB) <- node_names

# create the igraph object for the full network
wadj <- graph_from_adjacency_matrix(
  AB,
  mode = "directed",
  weighted = TRUE,
  diag = FALSE
)

# set vertex attributes for plants and herbivores
V(wadj)$color <- c(rep("grey50", nrow(alpha)), rep("white", 6))
V(wadj)$label.color <- c(rep("white", nrow(alpha)), rep("grey50", 6))
V(wadj)$frame.color <- c(rep("grey50", nrow(alpha)), rep("grey50", 6))

# assign edge curvature based on the number of outgoing edges from each node
cedge <- rep(0, length(E(wadj)$weight))
for (i in 1:vcount(wadj)) {
  edges_from_node <- as.numeric(E(wadj)[.from(i)])
  if (length(edges_from_node) > 0) {
    half <- round(floor(length(edges_from_node) / 2))
    cedge[edges_from_node[1:half]] <- -0.15
    cedge[edges_from_node[(half + 1):length(edges_from_node)]] <- 0.15
  }
}
E(wadj)$curv <- cedge

# create separate graphs for positive and negative interactions
wadj_pos <- delete_edges(wadj, E(wadj)[weight < 0])
wadj_neg <- delete_edges(wadj, E(wadj)[weight > 0])

# set edge attributes for positive and negative interaction graphs
E(wadj_pos)$color <- adjustcolor("#5e9ca0", alpha.f = 0.85)
E(wadj_pos)$lty <- 1
E(wadj_pos)$weight <- abs(E(wadj_pos)$weight)

E(wadj_neg)$color <- adjustcolor("#d4847b", alpha.f = 0.85)
E(wadj_neg)$lty <- 1
E(wadj_neg)$weight <- abs(E(wadj_neg)$weight)

# create a layout for the full network
set.seed(100)
colnames(alpha) <- colnames(AB)[1:nrow(alpha)]
wadj_alpha <- graph_from_adjacency_matrix(
  as.matrix(alpha),
  mode = "directed",
  weighted = TRUE,
  diag = FALSE
)
E(wadj_alpha)$weight <- abs(E(wadj_alpha)$weight)
lyo_alpha <- layout.sphere(wadj_alpha)
lyo <- layout.sphere(wadj)
lyo[1:nrow(alpha), ] <- lyo_alpha

# manually adjust the layout to spread out the nodes and reduce overlap
lyo[nrow(alpha) + 1, ] <- c(-0.90, 1.1, 1)
lyo[nrow(alpha) + 2, ] <- c(-0.65, 1.5, 1)
lyo[nrow(alpha) + 3, ] <- c(-0.25, 1.7, 1)
lyo[nrow(alpha) + 4, ] <- c(0.25, 1.7, 1)
lyo[nrow(alpha) + 5, ] <- c(0.65, 1.5, 1)
lyo[nrow(alpha) + 6, ] <- c(0.90, 1.1, 1)

adjustments <- list(
  "Am" = c(0.10, 0.10),
  "Vp" = c(-0.10, -0.10),
  "Pt" = c(-0.15, 0.00),
  "Pp" = c(-0.10, 0.10),
  "Tf" = c(-0.10, 0.10),
  "Pa" = c(-0.05, -0.10),
  "Pl" = c(-0.05, 0.00),
  "Dg" = c(0.00, -0.10),
  "Ra" = c(-0.10, 0.00),
  "Dc" = c(0.10, -0.10),
  "Er" = c(0.15, 0.00),
  "Em" = c(0.05, -0.10),
  "Fa" = c(0.05, -0.05),
  "Pe" = c(0.05, 0.05),
  "So" = c(0.05, 0.00),
  "Be" = c(0.075, 0.00),
  "Vo" = c(-0.20, -0.225),
  "Ae" = c(0.10, 0.15),
  "Fr" = c(0.05, 0.00),
  "Cj" = c(0.05, 0.05),
  "Ca" = c(-0.025, -0.025)
)
for (node in names(adjustments)) {
  idx <- which(colnames(AB) == node)
  lyo[idx, 1:2] <- lyo[idx, 1:2] + adjustments[[node]]
}

# plot the networks and save to a temporary file
temp_network_file <- tempfile(fileext = ".png")
png(temp_network_file, width = 20, height = 10, units = "in", res = 500)
par(mfrow = c(1, 2), mar = c(3, 1, 3, 1))

# plot the positive interaction network
plot(
  wadj_pos,
  layout = lyo,
  edge.color = E(wadj_pos)$color,
  edge.curved = E(wadj_pos)$curv,
  edge.width = (E(wadj_pos)$weight * 1.5 + 1.5),
  edge.arrow.size = 0.75,
  vertex.shape = "circle",
  vertex.size = 13,
  vertex.label.family = "sans",
  vertex.label.color = V(wadj_pos)$label.color,
  vertex.label.font = 2,
  vertex.label.cex = 1.75,
  vertex.color = V(wadj_pos)$color,
  vertex.frame.color = V(wadj_pos)$frame.color
)

# add a custom legend to the centre of both networks
legend(
  x = par("usr")[2] + 0.05,
  y = par("usr")[3] + 0.2,
  xjust = 0.5,
  yjust = 0.6,
  xpd = NA,
  legend = c(
    "Positive interaction",
    "Negative interaction",
    "Plants",
    "Grasshoppers"
  ),
  col = c("#5e9ca0", "#d4847b", "grey50", "grey50"),
  pt.bg = c(NA, NA, "grey50", "white"),
  lwd = c(4, 4, NA, NA),
  pch = c(NA, NA, 21, 21),
  pt.cex = c(NA, NA, 3, 3),
  bty = "n",
  cex = 1.6,
  text.font = 2
)

# plot the negative interaction network
plot(
  wadj_neg,
  layout = lyo,
  edge.color = E(wadj_neg)$color,
  edge.curved = E(wadj_neg)$curv,
  edge.width = (E(wadj_neg)$weight * 1.5 + 1.5),
  edge.arrow.size = 0.75,
  vertex.shape = "circle",
  vertex.size = 13,
  vertex.label.family = "sans",
  vertex.label.color = V(wadj_neg)$label.color,
  vertex.label.font = 2,
  vertex.label.cex = 1.75,
  vertex.color = V(wadj_neg)$color,
  vertex.frame.color = V(wadj_neg)$frame.color
)
invisible(dev.off())


# interaction layers histograms ------------------------------------------

# combine the alpha, gamma, and beta matrices
beta <- Reduce('+', beta_gp) + Reduce('+', beta_pp)

# helper function to extract non-zero coefficients and label them by type
extract_coefs <- function(mat, name) {
  vals <- as.vector(mat)
  vals <- vals[vals != 0 & !is.na(vals)]
  data.frame(coef = vals, type = name)
}

# helper function to calculate connectance (proportion of non-zero interactions)
calc_connectance <- function(mat) {
  non_zeros <- sum(mat != 0, na.rm = TRUE)
  total <- length(as.vector(mat))
  round(non_zeros / total, 2)
}

# extract coefficients and calculate connectance for each layer
cdist <- bind_rows(
  extract_coefs(alpha, "alpha"),
  extract_coefs(gamma, "gamma"),
  extract_coefs(beta, "beta")
)
cdist$type <- factor(cdist$type, levels = c("alpha", "gamma", "beta"))

# create a data frame for connectance labels to be added to the plot
con_text <- data.frame(
  type = factor(
    c("alpha", "gamma", "beta"),
    levels = c("alpha", "gamma", "beta")
  ),
  connectance = c(
    paste("Cn =", calc_connectance(alpha)),
    paste("Cn =", calc_connectance(gamma)),
    paste("Cn =", calc_connectance(beta))
  )
)

# define colors and labels for the layers to be used in the plot
layer_colors <- c("alpha" = "grey20", "gamma" = "grey50", "beta" = "grey75")
layer_labels <- c(
  "alpha" = expression(alpha),
  "gamma" = expression(gamma),
  "beta" = expression(beta)
)

# create the histogram and density plot for the interaction coefficients
plot_histograms <- ggplot(cdist, aes(x = coef)) +
  geom_histogram(
    aes(y = after_stat(density), fill = type, color = type),
    binwidth = 0.005,
    alpha = 0.75
  ) +
  geom_density(
    aes(fill = type, color = type),
    outline.type = "full",
    alpha = 0.5
  ) +
  geom_vline(
    xintercept = 0,
    linetype = "dotted",
    color = "black",
    size = 1.25
  ) +
  geom_text(
    data = con_text,
    mapping = aes(x = -0.48, y = Inf, label = connectance),
    vjust = 1.5,
    hjust = 0,
    size = 7.5,
    fontface = "bold",
    family = "sans"
  ) +
  facet_grid(rows = vars(type), scales = "free_y") +
  scale_fill_manual(
    name = "Network layers:",
    values = layer_colors,
    labels = layer_labels
  ) +
  scale_color_manual(
    name = "Network layers:",
    values = layer_colors,
    labels = layer_labels
  ) +
  labs(x = "Interaction coefficients", y = "Density") +
  coord_cartesian(xlim = c(-0.5, 0.5)) +
  theme_classic() +
  theme(
    legend.direction = "horizontal",
    legend.position = c(0.85, 0.90),
    # Scaled up legend text a tiny bit more
    legend.title = element_text(face = "bold", size = 20),
    legend.text = element_text(size = 20, face = "bold"),
    # Scaled up axis titles and text to match the network's visual weight
    axis.title = element_text(
      size = 20,
      face = "bold",
      margin = margin(t = 12, r = 12)
    ),
    axis.text = element_text(size = 20, color = "black"),
    strip.background = element_blank(),
    strip.text.y = element_blank()
  )


# combine plots ----------------------------------------------------------

# read the saved network image and convert it to a ggplot object
p_networks_clean <- ggdraw() + draw_image(temp_network_file)

# create a spacer plot to add vertical space between the network and histogram
spacer <- ggplot() + theme_void()
bottom_row <- plot_grid(
  spacer,
  plot_histograms,
  spacer,
  rel_widths = c(0.05, 0.9, 0.05),
  nrow = 1
)

# combine the network and histogram plots
final_figure <- plot_grid(
  p_networks_clean,
  spacer, # gap between top and bottom halves
  bottom_row,
  ncol = 1,
  rel_heights = c(1, 0.05, 0.8), # allocated 5% height to the vertical gap
  labels = c("a", "", "b"), # empty string for the spacer's label
  label_size = 32,
  label_x = 0.02,
  label_y = 1,
  label_fontface = "bold"
)

# export -----------------------------------------------------------------

# save the final combined figure as a high-resolution JPEG
ggsave(
  filename = "suppl-mat/suppl-figures/fig-networks-histograms.jpeg",
  plot = final_figure,
  device = "jpeg",
  dpi = 640,
  width = 16,
  height = 16,
  bg = "white"
)

# clean up temp file
unlink(temp_network_file)

# info
cat(
  "Network and histogram plots have been successfully created and saved.\n"
)
