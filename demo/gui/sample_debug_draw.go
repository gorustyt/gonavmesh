package gui

type SampleTool interface {
	Type() int
	init(sample *Sample)
	reset()
	handleMenu()
	handleClick(s []float64, p []float64, shift bool)
	handleRender()
	handleRenderOverlay(proj, model *float64, view []int)
	handleToggle()
	handleStep()
	handleUpdate(dt float64)
}

type SampleToolState interface {
	init(sample *Sample)
	reset()
	handleRender()
	handleRenderOverlay(proj, model *float64, view *int)
	handleUpdate(dt float64)
}
