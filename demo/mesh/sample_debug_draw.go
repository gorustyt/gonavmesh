package mesh

type ISample interface {
	handleMeshChanged(geom *InputGeom)
	setTool(tool SampleTool)
	HandleRender()
}
type SampleTool interface {
	Type() int
	init(sample *Sample)
	reset()
	handleClick(s []float32, p []float32, shift bool)
	handleRender()
	handleRenderOverlay(proj, model []float32, view []int)
	handleToggle()
	handleStep()
	handleUpdate(dt float32)
}

type SampleToolState interface {
	init(sample *Sample)
	reset()
	handleRender()
	handleRenderOverlay(proj, model []float32, view []int)
	handleUpdate(dt float32)
}
