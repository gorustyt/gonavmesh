package debug_utils

type Colorb [4]uint8

func (c Colorb) R() uint8 {
	return c[0]
}

func (c Colorb) G() uint8 {
	return c[1]
}
func (c Colorb) Int() uint32 {
	return uint32(c.R()) | (uint32(c.G()) << 8) | (uint32(c.B()) << 16) | (uint32(c.A()) << 24)
}

func (c *Colorb) FromInt(col uint32) {
	c[0] = uint8(col & 0xff)
	c[1] = uint8((col >> 8) & 0xff)
	c[2] = uint8((col >> 16) & 0xff)
	c[3] = uint8((col >> 24) & 0xff)
}

func (c Colorb) B() uint8 {
	return c[2]
}

func (c Colorb) A() uint8 {
	return c[3]
}
