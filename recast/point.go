package recast

type Point struct {
	X, Y, Z float64
}

func (p *Point) Copy() *Point {
	return &Point{X: p.X, Y: p.Y, Z: p.Z}
}
