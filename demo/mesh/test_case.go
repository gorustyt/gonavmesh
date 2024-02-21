package mesh

import (
	"bufio"
	"fmt"
	"github.com/go-gl/gl/v4.1-core/gl"
	"github.com/gorustyt/gonavmesh/common"
	"github.com/gorustyt/gonavmesh/recast"
	"io"
	"log"
	"os"
	"strings"
	"time"
)

type TestType int

const (
	TEST_PATHFIND TestType = iota
	TEST_RAYCAST
)

type Test struct {
	Type         TestType
	spos         []float64
	epos         []float64
	nspos        []float64
	nepos        []float64
	radius       float64
	includeFlags int
	excludeFlags int
	expand       bool

	straight  []float64
	nstraight int
	polys     []detour.DtPolyRef
	npolys    int

	findNearestPolyTime  int64
	findPathTime         int64
	findStraightPathTime int64

	next *Test
}

func newTest() *Test {
	return &Test{
		spos:  make([]float64, 3),
		epos:  make([]float64, 3),
		nspos: make([]float64, 3),
		nepos: make([]float64, 3),
	}
}

type TestCase struct {
	m_sampleName   string
	m_geomFileName string
	m_tests        *Test
	gs             *guiState
}

func newTestCase(gs *guiState) *TestCase {
	return &TestCase{gs: gs}
}

func (t *TestCase) parseRow(ss []string) {
	switch ss[0] {
	case "s":
		// Sample name.
		t.m_sampleName = ss[1]
	case "f":
		// File name.
		t.m_geomFileName = ss[1]
	case "pf":
		// Pathfind test.
		test := newTest()
		test.Type = TEST_PATHFIND
		test.expand = false
		test.next = t.m_tests
		t.m_tests = test
		_, err := fmt.Sscanf(strings.Join(ss[2:], " "), "%f %f %f %f %f %f %hx %hx",
			&test.spos[0], &test.spos[1], &test.spos[2],
			&test.epos[0], &test.epos[1], &test.epos[2],
			&test.includeFlags, &test.excludeFlags)
		if err != nil {
			panic(err)
		}
	case "rc":
		// Pathfind test.
		test := newTest()
		test.Type = TEST_RAYCAST
		test.expand = false
		t.m_tests.next = t.m_tests
		t.m_tests = test
		_, err := fmt.Sscanf(strings.Join(ss[2:], ""), "%f %f %f %f %f %f %hx %hx",
			&test.spos[0], &test.spos[1], &test.spos[2],
			&test.epos[0], &test.epos[1], &test.epos[2],
			&test.includeFlags, &test.excludeFlags)
		if err != nil {
			panic(err)
		}
	}
}
func (t *TestCase) load(path string) bool {
	f, err := os.Open(path)
	if err != nil {
		log.Println(err)
		panic("")
	}
	defer f.Close()
	reader := bufio.NewReader(f)
	line := ""
	for {
		lines, prefix, err := reader.ReadLine()
		if err == io.EOF {
			break
		}
		if err != nil {
			log.Println(err)
			panic("")
		}
		if prefix {
			line += string(lines)
		} else {
			t.parseRow(strings.Split(line, " "))
		}
	}

	return true
}

func (t *TestCase) getSampleName() string   { return t.m_sampleName }
func (t *TestCase) getGeomFileName() string { return t.m_geomFileName }
func (t *TestCase) resetTimes() {
	for iter := t.m_tests; iter != nil; iter = iter.next {
		iter.findNearestPolyTime = 0
		iter.findPathTime = 0
		iter.findStraightPathTime = 0
	}
}
func (t *TestCase) doTests(navmesh recast.IDtNavMesh, navquery recast.NavMeshQuery) {
	if navmesh == nil || navquery == nil {
		return
	}

	t.resetTimes()

	const MAX_POLYS = 256
	polys := make([]detour.DtPolyRef, MAX_POLYS)
	straight := make([]float64, MAX_POLYS*3)
	polyPickExt := []float64{2, 4, 2}

	for iter := t.m_tests; iter != nil; iter = iter.next {
		iter.polys = nil
		iter.npolys = 0
		iter.straight = nil
		iter.nstraight = 0

		var filter detour.DtQueryFilter
		filter.SetIncludeFlags(iter.includeFlags)
		filter.SetExcludeFlags(iter.excludeFlags)

		// Find start points
		findNearestPolyStart := time.Now()

		var startRef, endRef detour.DtPolyRef
		startRef, _ = navquery.FindNearestPoly(iter.spos[:], polyPickExt, &filter, iter.nspos[:])
		endRef, _ = navquery.FindNearestPoly(iter.epos[:], polyPickExt, &filter, iter.nepos[:])

		findNearestPolyEnd := time.Since(findNearestPolyStart)
		iter.findNearestPolyTime += findNearestPolyEnd.Milliseconds()

		if startRef == 0 || endRef == 0 {
			continue
		}

		if iter.Type == TEST_PATHFIND {
			// Find path
			findPathStart := time.Now()

			iter.npolys, _ = navquery.FindPath(startRef, endRef, iter.spos[:], iter.epos[:], &filter, polys, MAX_POLYS)

			findPathEnd := time.Since(findPathStart)
			iter.findPathTime += findPathEnd.Milliseconds()

			// Find straight path
			if iter.npolys > 0 {
				findStraightPathStart := time.Now()

				iter.nstraight, _ = navquery.FindStraightPath(iter.spos[:], iter.epos[:], polys, iter.npolys,
					straight, []int{}, []detour.DtPolyRef{}, MAX_POLYS)
				iter.findStraightPathTime += time.Since(findStraightPathStart).Milliseconds()
			}

			// Copy results
			if iter.npolys > 0 {
				iter.polys = make([]detour.DtPolyRef, iter.npolys)
				copy(iter.polys, polys[:iter.npolys])
			}
			if iter.nstraight > 0 {
				iter.straight = make([]float64, iter.nstraight*3)
				copy(iter.straight, straight[:3*iter.nstraight])
			}
		} else if iter.Type == TEST_RAYCAST {
			tt := 0.0

			hitNormal := make([]float64, 2*3)
			hitPos := make([]float64, 2*3)
			iter.straight = make([]float64, 2*3)
			iter.nstraight = 2

			iter.straight[0] = iter.spos[0]
			iter.straight[1] = iter.spos[1]
			iter.straight[2] = iter.spos[2]

			findPathStart := time.Now()

			navquery.Raycast(startRef, iter.spos[:], iter.epos[:], &filter, &tt, hitNormal, polys, &iter.npolys, MAX_POLYS)
			iter.findPathTime += time.Since(findPathStart).Milliseconds()

			if tt > 1 {
				// No hit
				copy(hitPos, iter.epos[:])
			} else {
				// Hit
				common.Vlerp(hitPos, iter.spos[:], iter.epos[:], tt)
			}
			// Adjust height.
			if iter.npolys > 0 {

				h, _ := navquery.GetPolyHeight(polys[iter.npolys-1], hitPos)
				hitPos[1] = h
			}
			copy(iter.straight[3:], hitPos)

			if iter.npolys > 0 {
				iter.polys = make([]detour.DtPolyRef, iter.npolys)
				copy(iter.polys, polys[:iter.npolys])
			}
		}
	}

	log.Printf("Test Results:\n")
	n := 0
	for iter := t.m_tests; iter != nil; iter = iter.next {
		total := iter.findNearestPolyTime + iter.findPathTime + iter.findStraightPathTime
		log.Printf(" - Path %02d:     %.4f ms\n", n, float64(total)/1000.0)
		log.Printf("    - poly:     %.4f ms\n", float64(iter.findNearestPolyTime)/1000.0)
		log.Printf("    - path:     %.4f ms\n", float64(iter.findPathTime)/1000.0)
		log.Printf("    - straight: %.4f ms\n", float64(iter.findStraightPathTime)/1000.0)
		n++
	}
}

func (t *TestCase) handleRender() {
	gl.LineWidth(2.0)
	glBegin(GL_LINES)
	for iter := t.m_tests; iter != nil; iter = iter.next {
		dir := make([]float64, 3)
		common.Vsub(dir, iter.epos[:], iter.spos[:])
		common.Vnormalize(dir)
		glColor4ub(128, 25, 0, 192)
		glVertex3f(iter.spos[0], iter.spos[1]-0.3, iter.spos[2])
		glVertex3f(iter.spos[0], iter.spos[1]+0.3, iter.spos[2])
		glVertex3f(iter.spos[0], iter.spos[1]+0.3, iter.spos[2])
		glVertex3f(iter.spos[0]+dir[0]*0.3, iter.spos[1]+0.3+dir[1]*0.3, iter.spos[2]+dir[2]*0.3)
		glColor4ub(51, 102, 0, 129)
		glVertex3f(iter.epos[0], iter.epos[1]-0.3, iter.epos[2])
		glVertex3f(iter.epos[0], iter.epos[1]+0.3, iter.epos[2])

		if iter.expand {
			s := 0.1
			glColor4ub(255, 32, 0, 128)
			glVertex3f(iter.spos[0]-s, iter.spos[1], iter.spos[2])
			glVertex3f(iter.spos[0]+s, iter.spos[1], iter.spos[2])
			glVertex3f(iter.spos[0], iter.spos[1], iter.spos[2]-s)
			glVertex3f(iter.spos[0], iter.spos[1], iter.spos[2]+s)
			glColor4ub(255, 192, 0, 255)
			glVertex3f(iter.nspos[0]-s, iter.nspos[1], iter.nspos[2])
			glVertex3f(iter.nspos[0]+s, iter.nspos[1], iter.nspos[2])
			glVertex3f(iter.nspos[0], iter.nspos[1], iter.nspos[2]-s)
			glVertex3f(iter.nspos[0], iter.nspos[1], iter.nspos[2]+s)

			glColor4ub(255, 32, 0, 128)
			glVertex3f(iter.epos[0]-s, iter.epos[1], iter.epos[2])
			glVertex3f(iter.epos[0]+s, iter.epos[1], iter.epos[2])
			glVertex3f(iter.epos[0], iter.epos[1], iter.epos[2]-s)
			glVertex3f(iter.epos[0], iter.epos[1], iter.epos[2]+s)
			glColor4ub(255, 192, 0, 255)
			glVertex3f(iter.nepos[0]-s, iter.nepos[1], iter.nepos[2])
			glVertex3f(iter.nepos[0]+s, iter.nepos[1], iter.nepos[2])
			glVertex3f(iter.nepos[0], iter.nepos[1], iter.nepos[2]-s)
			glVertex3f(iter.nepos[0], iter.nepos[1], iter.nepos[2]+s)
		}

		if iter.expand {
			glColor4ub(255, 192, 0, 255)
		} else {
			glColor4ub(0, 0, 0, 64)
		}

		for i := 0; i < iter.nstraight-1; i++ {
			glVertex3f(iter.straight[i*3+0], iter.straight[i*3+1]+0.3, iter.straight[i*3+2])
			glVertex3f(iter.straight[(i+1)*3+0], iter.straight[(i+1)*3+1]+0.3, iter.straight[(i+1)*3+2])
		}
	}
	gl.End()
	gl.LineWidth(1.0)
}
func (t *TestCase) handleRenderOverlay(proj, model []float64, view []int) bool {
	var text, subtext string
	n := 0

	LABEL_DIST := 1.0

	for iter := t.m_tests; iter != nil; iter = iter.next {
		pt := make([]float64, 3)
		dir := make([]float64, 3)
		if iter.nstraight > 0 {
			copy(pt, iter.straight[3:])
			if common.Vdist(pt, iter.spos[:]) > LABEL_DIST {
				common.Vsub(dir, pt, iter.spos[:])
				common.Vnormalize(dir)
				common.Vmad(pt, iter.spos[:], dir, LABEL_DIST)
			}
			pt[1] += 0.5
		} else {
			common.Vsub(dir, iter.epos[:], iter.spos[:])
			common.Vnormalize(dir)
			common.Vmad(pt, iter.spos[:], dir, LABEL_DIST)
			pt[1] += 0.5
		}
		res := common.GluProject([]float64{pt[0], pt[1], pt[2]},
			model, proj, view)
		if len(res) > 0 {
			x, y := int(res[0]), int(res[1])
			text = fmt.Sprintf("Path %d\n", n)
			col := imguiRGBA(0, 0, 0, 128)
			if iter.expand {
				col = imguiRGBA(255, 192, 0, 220)
			}

			t.gs.imguiDrawText(x, (y - 25), IMGUI_ALIGN_CENTER, text, col)
		}
		n++
	}

	resScroll := 0
	mouseOverMenu := t.gs.imguiBeginScrollArea("Test Results", 10, view[3]-10-350, 200, 350, &resScroll)
	//		mouseOverMenu = true;

	n = 0
	for iter := t.m_tests; iter != nil; iter = iter.next {
		total := iter.findNearestPolyTime + iter.findPathTime + iter.findStraightPathTime
		subtext = fmt.Sprintf("%.4f ms", float64(total)/1000.0)
		text = fmt.Sprintf("Path %d", n)

		if t.gs.imguiCollapse(text, subtext, iter.expand) {
			iter.expand = !iter.expand
		}

		if iter.expand {
			text = fmt.Sprintf("Poly: %.4f ms", float64(iter.findNearestPolyTime)/1000.0)
			t.gs.imguiValue(text)

			text = fmt.Sprintf("Path: %.4f ms", float64(iter.findPathTime)/1000.0)
			t.gs.imguiValue(text)

			text = fmt.Sprintf("Straight: %.4f ms", float64(iter.findStraightPathTime)/1000.0)
			t.gs.imguiValue(text)

			t.gs.imguiSeparator()
		}

		n++
	}

	t.gs.imguiEndScrollArea()

	return mouseOverMenu
}
