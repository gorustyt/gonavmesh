package mesh

import (
	"github.com/go-gl/mathgl/mgl32"
	"github.com/gorustyt/fyne/v2/canvas3d/context"
	"github.com/gorustyt/fyne/v2/canvas3d/context/enum"
	"github.com/gorustyt/gonavmesh/debug_utils"
	"github.com/gorustyt/gonavmesh/demo/config"
)

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

var (
	_ debug_utils.DuDebugDraw = (*DebugDrawGL)(nil)
)

type DebugDrawGL struct {
	ctx context.Context
}

func (g *DebugDrawGL) AreaToCol(area int) debug_utils.Colorb {
	switch area {
	// Ground (0) : light blue
	case config.SAMPLE_POLYAREA_GROUND:
		return debug_utils.DuRGBA(0, 192, 255, 255)
	// Water : blue
	case config.SAMPLE_POLYAREA_WATER:
		return debug_utils.DuRGBA(0, 0, 255, 255)
	// Road : brown
	case config.SAMPLE_POLYAREA_ROAD:
		return debug_utils.DuRGBA(50, 20, 12, 255)
	// Door : cyan
	case config.SAMPLE_POLYAREA_DOOR:
		return debug_utils.DuRGBA(0, 255, 255, 255)
	// Grass : green
	case config.SAMPLE_POLYAREA_GRASS:
		return debug_utils.DuRGBA(0, 255, 0, 255)
	// Jump : yellow
	case config.SAMPLE_POLYAREA_JUMP:
		return debug_utils.DuRGBA(255, 255, 0, 255)
	// Unexpected : red
	default:
		return debug_utils.DuRGBA(255, 0, 0, 255)
	}
}

func NewDebugDrawGL(ctx context.Context) debug_utils.DuDebugDraw {
	return &DebugDrawGL{
		ctx: ctx,
	}
}

var (
	g_tex GLCheckerTexture
)

func (g *DebugDrawGL) DepthMask(state bool) {
	if state {
		g.ctx.ExtDepthMask(enum.GlTrue != 0)
	} else {
		g.ctx.ExtDepthMask(enum.GlFalse != 0)
	}
}

func (g *DebugDrawGL) Texture(state bool) {
	if state {
		g.ctx.Enable(enum.Texture2D)
		g_tex.Bind()
	} else {
		g.ctx.Disable(enum.Texture2D)
	}
}

func (g *DebugDrawGL) Begin(prim debug_utils.DuDebugDrawPrimitives, sizes ...float32) {
	size := float32(1.)
	if len(sizes) > 0 {
		size = sizes[0]
	}
	switch prim {
	case debug_utils.DU_DRAW_POINTS:
		g.ctx.ExtPointSize(size)
		g.ctx.ExtBegin(enum.Points)
		break
	case debug_utils.DU_DRAW_LINES:
		g.ctx.ExtLineWidth(size)
		g.ctx.ExtBegin(enum.Lines)
		break
	case debug_utils.DU_DRAW_TRIS:
		g.ctx.ExtBegin(enum.Triangles)
		break
	case debug_utils.DU_DRAW_QUADS:
		g.ctx.ExtBegin(enum.Quads)
		break
	}
}

func (g *DebugDrawGL) Vertex(pos []float32, color debug_utils.Colorb) {
	g.ctx.ExtColor4ubv(color[:])
	g.ctx.ExtVertex3fv(mgl32.Vec3{pos[0], pos[1], pos[2]})
}

func (g *DebugDrawGL) Vertex1(x, y, z float32, color debug_utils.Colorb) {
	g.ctx.ExtColor4ubv(color[:])
	g.ctx.ExtVertex3f(x, y, z)
}

func (g *DebugDrawGL) Vertex2(pos []float32, color debug_utils.Colorb, uv []float32) {
	g.ctx.ExtColor4ubv(color[:])
	g.ctx.ExtTexCoord2fv(mgl32.Vec2{uv[0], uv[1]})
	g.ctx.ExtVertex3fv(mgl32.Vec3{pos[0], pos[1], pos[2]})
}

func (g *DebugDrawGL) Vertex3(x, y, z float32, color debug_utils.Colorb, u, v float32) {
	g.ctx.ExtColor4ubv(color[:])
	g.ctx.ExtTexCoord2fv(mgl32.Vec2{u, v})
	g.ctx.ExtVertex3f(x, y, z)
}

func (g *DebugDrawGL) End() {
	g.ctx.ExtEnd()
	g.ctx.ExtLineWidth(1.0)
	g.ctx.ExtPointSize(1.0)
}

type GLCheckerTexture struct {
	texId uint32
	ctx   context.Context
}

func (c *GLCheckerTexture) Bind() {
	if c.texId == 0 {
		// Create checker pattern.
		col0 := []uint8{215, 215, 215, 255}
		col1 := []uint8{255, 255, 255, 255}
		const TSIZE = 64
		var data [TSIZE * TSIZE][]uint8
		c.ctx.GenTextures(1, &c.texId)
		c.ctx.BindTexture(enum.Texture2D, context.Texture(c.texId))
		level := 0
		size := TSIZE
		for size > 0 {
			for y := 0; y < size; y++ {
				for x := 0; x < size; x++ {
					col := col1
					if x == 0 || y == 0 {
						col = col0
					}
					data[x+y*size] = col
				}
			}
			c.ctx.TexImage2D(uint32(enum.Texture2D), level, size, size, enum.GlRgba, enum.GlUnsigedBytes, data[0]) //TODO verify
			size /= 2
			level++
		}
		c.ctx.TexParameteri(enum.Texture2D, enum.TextureMinFilter, enum.LinearMipMapNearest)
		c.ctx.TexParameteri(enum.Texture2D, enum.TextureMagFilter, enum.Linear)
	} else {
		c.ctx.BindTexture(enum.Texture2D, context.Texture(c.texId))
	}
}
