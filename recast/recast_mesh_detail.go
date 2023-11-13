package recast

const (
	RC_UNSET_HEIGHT = 0xffff
)

type rcHeightPatch struct {
	data                      []int
	xmin, ymin, width, height int
}
