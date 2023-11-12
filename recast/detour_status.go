package recast

type dtStatus int

const (
	// High level status.
	DT_FAILURE     = 1 << 31 // Operation failed.
	DT_SUCCESS     = 1 << 30 // Operation succeed.
	DT_IN_PROGRESS = 1 << 29 // Operation still in progress.

	// Detail information for status.
	DT_STATUS_DETAIL_MASK = 0x0ffffff
	DT_WRONG_MAGIC        = 1 << 0 // Input data is not recognized.
	DT_WRONG_VERSION      = 1 << 1 // Input data is in wrong version.
	DT_OUT_OF_MEMORY      = 1 << 2 // Operation ran out of memory.
	DT_INVALID_PARAM      = 1 << 3 // An input parameter was invalid.
	DT_BUFFER_TOO_SMALL   = 1 << 4 // Result buffer for the query was too small to store all results.
	DT_OUT_OF_NODES       = 1 << 5 // Query ran out of nodes during search.
	DT_PARTIAL_RESULT     = 1 << 6 // Query did not reach the end location, returning best guess.
	DT_ALREADY_OCCUPIED   = 1 << 7 // A tile has already been assigned to the given x,y coordinate
)

// Returns true of status is success.
func (status dtStatus) dtStatusSucceed() bool {
	return (status & DT_SUCCESS) != 0
}

// Returns true of status is failure.
func (status dtStatus) dtStatusFailed() bool {
	return (status & DT_FAILURE) != 0
}

// Returns true of status is in progress.
func (status dtStatus) dtStatusInProgress() bool {
	return (status & DT_IN_PROGRESS) != 0
}

// Returns true if specific detail is set.
func (status dtStatus) dtStatusDetail(detail int) bool {
	return (int(status) & detail) != 0
}
