package recast

// / Sets the neighbor connection data for the specified direction.
// / @param[in]		span			The span to update.
// / @param[in]		direction		The direction to set. [Limits: 0 <= value < 4]
// / @param[in]		neighborIndex	The index of the neighbor span.
func rcSetCon(span *RcCompactSpan, direction, neighborIndex int32) {
	shift := direction * 6
	con := span.Con
	span.Con = (con & ^(0x3f << shift)) | ((neighborIndex & 0x3f) << shift)
}
