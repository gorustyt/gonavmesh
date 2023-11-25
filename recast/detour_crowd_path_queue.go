package recast

import "gonavamesh/common"

type dtPathQueueRef int

const DT_PATHQ_INVALID = 0

type PathQuery struct {
	ref dtPathQueueRef
	/// Path find start and end location.
	startPos, endPos [3]float64
	startRef, endRef DtPolyRef
	/// Result.
	path  []DtPolyRef
	npath int
	/// State.
	status    DtStatus
	keepAlive int
	filter    *DtQueryFilter ///< TODO: This is potentially dangerous!
}

type dtPathQueue struct {
	MAX_QUEUE     int
	m_queue       []*PathQuery
	m_nextHandle  dtPathQueueRef
	m_maxPathSize int
	m_queueHead   int
	m_navquery    NavMeshQuery
}

func newDtPathQueue() *dtPathQueue {
	d := &dtPathQueue{
		MAX_QUEUE: 8,
	}
	d.m_queue = make([]*PathQuery, d.MAX_QUEUE)
	return d
}

func (d *dtPathQueue) purge() {
	d.m_navquery = nil
	for i := 0; i < d.MAX_QUEUE; i++ {
		d.m_queue[i].path = nil
	}
}

func (d *dtPathQueue) init(maxPathSize int, maxSearchNodeCount int, nav IDtNavMesh) bool {
	d.purge()
	d.m_navquery = NewDtNavMeshQuery(nav, maxSearchNodeCount)
	d.m_maxPathSize = maxPathSize
	for i := 0; i < d.MAX_QUEUE; i++ {
		d.m_queue[i].ref = DT_PATHQ_INVALID
		d.m_queue[i].path = make([]DtPolyRef, d.m_maxPathSize)
	}

	d.m_queueHead = 0

	return true
}

func (d *dtPathQueue) update(maxIters int) {
	MAX_KEEP_ALIVE := 2 // in update ticks.

	// Update path request until there is nothing to update
	// or upto maxIters pathfinder iterations has been consumed.
	iterCount := maxIters

	for i := 0; i < d.MAX_QUEUE; i++ {
		q := d.m_queue[d.m_queueHead%d.MAX_QUEUE]

		// Skip inactive requests.
		if q.ref == DT_PATHQ_INVALID {
			d.m_queueHead++
			continue
		}

		// Handle completed request.
		if q.status.DtStatusSucceed() || q.status.DtStatusFailed() {
			// If the path result has not been read in few frames, free the slot.
			q.keepAlive++
			if q.keepAlive > MAX_KEEP_ALIVE {
				q.ref = DT_PATHQ_INVALID
				q.status = 0
			}

			d.m_queueHead++
			continue
		}

		// Handle query start.
		if q.status == 0 {
			q.status = d.m_navquery.InitSlicedFindPath(q.startRef, q.endRef, q.startPos[:], q.endPos[:], q.filter, 0)
		}
		// Handle query in progress.
		if q.status.DtStatusInProgress() {
			iters := 0
			iters, q.status = d.m_navquery.UpdateSlicedFindPath(iterCount)
			iterCount -= iters
		}
		if q.status.DtStatusSucceed() {
			q.npath, q.status = d.m_navquery.FinalizeSlicedFindPath(q.path, d.m_maxPathSize)
		}

		if iterCount <= 0 {
			break
		}

		d.m_queueHead++
	}
}

func (d *dtPathQueue) request(startRef, endRef DtPolyRef,
	startPos, endPos []float64,
	filter *DtQueryFilter) dtPathQueueRef {
	// Find empty slot
	slot := -1
	for i := 0; i < d.MAX_QUEUE; i++ {
		if d.m_queue[i].ref == DT_PATHQ_INVALID {
			slot = i
			break
		}
	}
	// Could not find slot.
	if slot == -1 {
		return DT_PATHQ_INVALID
	}

	ref := d.m_nextHandle
	d.m_nextHandle++
	if d.m_nextHandle == DT_PATHQ_INVALID {
		d.m_nextHandle++
	}

	q := d.m_queue[slot]
	q.ref = ref
	copy(q.startPos[:], startPos)
	q.startRef = startRef
	copy(q.endPos[:], endPos)
	q.endRef = endRef

	q.status = 0
	q.npath = 0
	q.filter = filter
	q.keepAlive = 0

	return ref
}

func (d *dtPathQueue) getRequestStatus(ref dtPathQueueRef) DtStatus {
	for i := 0; i < d.MAX_QUEUE; i++ {
		if d.m_queue[i].ref == ref {
			return d.m_queue[i].status
		}

	}
	return DT_FAILURE
}

func (d *dtPathQueue) getPathResult(ref dtPathQueueRef, path []DtPolyRef, pathSize *int, maxPath int) DtStatus {
	for i := 0; i < d.MAX_QUEUE; i++ {
		if d.m_queue[i].ref == ref {
			q := d.m_queue[i]
			details := q.status & DT_STATUS_DETAIL_MASK
			// Free request for reuse.
			q.ref = DT_PATHQ_INVALID
			q.status = 0
			// Copy path
			n := common.Min(q.npath, maxPath)
			copy(path, q.path[:n])
			*pathSize = n
			return details | DT_SUCCESS
		}
	}
	return DT_FAILURE
}
