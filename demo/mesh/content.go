package mesh

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/gonavmesh/demo/config"
	"log/slog"
)

type Content struct {
	cfg    *config.Config
	sample ISample
	geom *InputGeom
}

func (c *Content) GetConfig() *config.Config {
	return c.cfg
}

func (c *Content) MinSize() fyne.Size {
	return fyne.Size{
		Width:  600,
		Height: 1080,
	}
}

func (c *Content) Move(position fyne.Position) {

}

func (c *Content) Position() fyne.Position {
	return fyne.Position{}
}

func (c *Content) Resize(size fyne.Size) {

}

func (c *Content) Size() fyne.Size {
	return fyne.Size{}
}

func (c *Content) Hide() {

}

func (c *Content) Visible() bool {
	return true
}

func (c *Content) Show() {

}

func (c *Content) Refresh() {

}

func (c *Content) InputMeshChange(){
	geom:=newInputGeom()
	if !geom.load(""){
		slog.Info("geom load error")
		return
	}
	c.geom=geom
	if c.sample!=nil{
		c.sample.handleMeshChanged(geom)
	}
}

func (c *Content) SampleChange(sample string){
	switch sample {
	case config.SampleSoloMesh:
		c.sample=newSampleSoloMesh(c)
	case config.SampleTileMesh:
		c.sample=newSampleTileMesh(c)
	case config.SampleTempObstacles:
		c.sample=newSampleTempObstacles(c)
	}
}
func NewContent() *Content {
	cfg := config.NewConfig()
	c := &Content{
		cfg: cfg}
	cfg.PropsConfig.OnInputMesh= c.InputMeshChange
	return c
}
