package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/data/binding"
	"github.com/gorustyt/fyne/v2/widget"
	"gonavamesh/demo/config"
)

type ToolsCrowds struct {
	c *fyne.Container
}

func NewToolsCrowds(ctx *Context) *ToolsCrowds {
	cfg := ctx.Config()
	c := &ToolsCrowds{}
	group := widget.NewRadioGroup([]string{
		config.DescCrowdTool_TOOLMODE_CREATE,
		config.DescCrowdTool_TOOLMODE_MOVE_TARGET,
		config.DescCrowdTool_TOOLMODE_SELECT,
		config.DescCrowdTool_TOOLMODE_TOGGLE_POLYS,
	}, func(s string) {
		cfg.ToolsConfig.ToolModel = s
	})

	s1 := widget.NewSliderWithData(0, 3, binding.BindFloat(&cfg.ToolsConfig.ObstacleAvoidanceType))
	s1.Step = 1
	s2 := widget.NewSliderWithData(0, 20, binding.BindFloat(&cfg.ToolsConfig.SeparationWeight))
	s2.Step = 0.01
	optionsContainer := container.NewVBox(
		widget.NewLabel("Options"),
		widget.NewCheckGroup([]string{
			config.ExpandOptionsOptimizeVisibility,
			config.ExpandOptionsOptimizeTopology,
			config.ExpandOptionsAnticipateTurns,
			config.ExpandOptionsObstacleAvoidance,
			config.ExpandOptionsSeparation,
		}, func(strings []string) {
			cfg.ToolsConfig.ExpandOptions = strings
		}),
		widget.NewLabel("Avoidance Quality"),
		s1,
		widget.NewLabel("Separation Weight"),
		s2,
	)
	selectDebugDraw := widget.NewCheckGroup([]string{
		config.ExpandSelectedDebugDrawShowCorners,
		config.ExpandSelectedDebugDrawShowCollisionSegs,
		config.ExpandSelectedDebugDrawShowPath,
		config.ExpandSelectedDebugDrawShowVO,
		config.ExpandSelectedDebugDrawShowPathOptimization,
		config.ExpandSelectedDebugDrawShowNeighbours,
	}, func(strings []string) {
		cfg.ToolsConfig.ExpandSelectedDebugDraw = strings
	})
	debugDrawGroup := widget.NewCheckGroup([]string{
		config.ExpandDebugDrawShowLabels,
		config.ExpandDebugDrawShowProxGrid,
		config.ExpandDebugDrawShowNodes,
		config.ExpandDebugDrawShowPerfGraph,
		config.ExpandDebugDrawShowDetailAll}, func(strings []string) {
		cfg.ToolsConfig.ExpandDebugDraw = strings
	})
	c.c = container.NewVBox(
		group,
		widget.NewSeparator(),
		optionsContainer,
		widget.NewSeparator(),
		widget.NewLabel("Selected Debug Draw"),
		selectDebugDraw,
		widget.NewSeparator(),
		widget.NewLabel("Debug Draw"),
		debugDrawGroup,
	)
	return c
}

func (t *ToolsCrowds) GetRenderObjs() (res []fyne.CanvasObject) {
	return append(res, t.c)
}
