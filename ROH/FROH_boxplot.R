#finding froh
nroh_sroh$FROH <- nroh_sroh$TOTAL/1152448282

#box plot
ggplot(nroh_sroh, aes(x = group, y = FROH, fill = group)) + 
  geom_boxplot() +
  scale_fill_manual(values = c("BBIP" = "#61bab8", "SBJO" = "#deef7f", "SJVE" = "#faa56f")) +
  scale_y_continuous(
    breaks = seq(0, 0.5, by = 0.05),
    expand = expansion(mult = c(0.05, 0.1))
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

