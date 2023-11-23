package gui

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"os"
	"path"
	"strconv"
	"strings"
	"testing"
)

func assertTrue(t *testing.T, value bool, msg string) {
	if !value {
		t.Errorf(msg)
	}
}

func Compare[T int | float64](t *testing.T, real []T, expect string) {
	ss := strings.Split(expect, ",")
	for ss[len(ss)-1] == "" {
		ss = ss[:len(ss)-1]
	}
	assertTrue(t, len(ss) == len(real), "")
	eps := 0.009
	for i, v := range ss {
		v1 := fmt.Sprintf("%v", real[i])
		if j := strings.Index(v, "."); j != -1 {
			f, err := strconv.ParseFloat(v, 64)
			if err != nil {
				panic(err)
			}
			d := math.Abs(f - float64(real[i]))
			assertTrue(t, d <= eps, "")
		} else {
			assertTrue(t, v1 == v, "")
		}

	}

}

func CompareData(t *testing.T, p string, m *mesh) {
	p = path.Join("./bin/TestLoadData", fmt.Sprintf("%v.txt", p))
	var m_normals, m_verts, m_tris string
	f, err := os.Open(p)
	if err != nil {
		panic(err)
	}
	defer f.Close()
	lineCount := 0
	reader := bufio.NewReader(f)
	for {

		line, isPrefix, err := reader.ReadLine()
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		switch lineCount {
		case 0:
			m_normals += string(line)
		case 1:
			m_verts += string(line)

		case 2:
			m_tris += string(line)

		}
		if isPrefix {
			continue
		}
		lineCount++
	}
	assertTrue(t, lineCount == 3, "")
	Compare(t, m.m_normals, m_normals)
	Compare(t, m.m_verts, m_verts)
	Compare(t, m.m_tris, m_tris)
}

func TestMesh(t *testing.T) {
	s := mewMesh()
	s.Load("./bin/Meshes/nav_test.obj")
	assertTrue(t, s.m_vertCount == 884, "")
	assertTrue(t, s.m_triCount == 1612, "")
	CompareData(t, "nav_test.obj", s)

	s = mewMesh()
	s.Load("./bin/Meshes/dungeon.obj")
	assertTrue(t, s.m_vertCount == 5101, "")
	assertTrue(t, s.m_triCount == 10133, "")
	CompareData(t, "dungeon.obj", s)

	s = mewMesh()
	s.Load("./bin/Meshes/undulating.obj")
	assertTrue(t, s.m_vertCount == 2800, "")
	assertTrue(t, s.m_triCount == 5202, "")
	CompareData(t, "undulating.obj", s)
}
