df = df_sel_ratio

set_default_plot_size(14cm, 14cm)
p = Gadfly.plot(
    df,
    x=:x,
    y=:sel_ratio,
    Geom.smooth(smoothing=0.8),
    Guide.xlabel("Air Temperature (°C)"),
    Guide.ylabel("selectivity_ratio")
)

draw(PDF("./output/$(mode)/plot_base_$(mode).pdf", 14cm, 14cm), p)
draw(PNG("./output/$(mode)/plot_base_$(mode).png", 14cm, 14cm), p)