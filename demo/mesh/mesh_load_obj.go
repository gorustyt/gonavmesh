package mesh

import (
	"bufio"
	"io"
	"math"
	"os"
	"path"
	"strconv"
	"strings"
)

// 负责加载对象数据
type rcMeshLoaderObj struct {
	m_filename  string
	m_scale     float32
	m_verts     []float32
	m_tris      []int
	m_normals   []float32
	m_vertCount int
	m_triCount  int
}

func newRcMeshLoaderObj() *rcMeshLoaderObj {
	return &rcMeshLoaderObj{}
}
func (m *rcMeshLoaderObj) getVerts() []float32   { return m.m_verts }
func (m *rcMeshLoaderObj) getNormals() []float32 { return m.m_normals }
func (m *rcMeshLoaderObj) getTris() []int        { return m.m_tris }
func (m *rcMeshLoaderObj) getVertCount() int     { return m.m_vertCount }
func (m *rcMeshLoaderObj) getTriCount() int      { return m.m_triCount }
func (m *rcMeshLoaderObj) getFileName() string   { return m.m_filename }
func mewMesh() *rcMeshLoaderObj {
	return &rcMeshLoaderObj{m_scale: 1}
}

func (m *rcMeshLoaderObj) Load(p string) error {
	f, err := os.Open(p)
	defer f.Close()
	if err != nil {
		return err
	}
	reader := bufio.NewReader(f)
	for {
		line, _, err := reader.ReadLine()
		row := strings.TrimSpace(string(line))
		if err == io.EOF {
			break
		}
		if err != nil {
			return err
		}
		if strings.HasSuffix(row, "#") {
			continue
		}
		m.parseRow(strings.Split(row, " "))
	}
	m.m_filename = path.Base(p)
	return nil
}
func (m *rcMeshLoaderObj) parseVertex(ss []string) {
	count := 0
	for count < len(ss) {
		x, err := strconv.ParseFloat(ss[count], 64)
		if err != nil {
			panic(err)
		}
		count++
		y, err := strconv.ParseFloat(ss[count], 64)
		if err != nil {
			panic(err)
		}
		count++
		z, err := strconv.ParseFloat(ss[count], 64)
		if err != nil {
			panic(err)
		}
		m.addVertex(x, y, z)
		count++
	}
}

func (m *rcMeshLoaderObj) parseFace(ss []string) {
	count := 0
	getV := func(v string) int {
		vi, err := strconv.Atoi(v)
		if err != nil {
			panic(err)
		}
		if vi < 0 {
			return vi + m.m_vertCount
		}
		return vi - 1
	}
	data := []int{}
	for count < len(ss) {
		vs := strings.Split(ss[count], "/")
		data = append(data, getV(vs[len(vs)-1]))
		if len(data) > 32 {
			break
		}
		count++
	}
	for i := 2; i < len(data); i++ {
		a := data[0]
		b := data[i-1]
		c := data[i]
		if a < 0 || a >= m.m_vertCount || b < 0 || b >= m.m_vertCount || c < 0 || c >= m.m_vertCount {
			continue
		}
		m.addTriangle(a, b, c)
	}
}
func (m *rcMeshLoaderObj) parseRow(ss []string) {

	switch ss[0] {
	case "v":
		if ss[1] != "n" && ss[1] != "t" {
			m.parseVertex(ss[1:])
		}
	case "f":
		m.parseFace(ss[1:])
	}
	// Calculate normals.
	m.m_normals = make([]float32, m.m_triCount*3)
	for i := 0; i < m.m_triCount*3; i += 3 {
		v0 := m.m_verts[m.m_tris[i]*3:]
		v1 := m.m_verts[m.m_tris[i+1]*3:]
		v2 := m.m_verts[m.m_tris[i+2]*3:]
		e0 := make([]float32, 3)
		e1 := make([]float32, 3)
		for j := 0; j < 3; j++ {
			e0[j] = v1[j] - v0[j]
			e1[j] = v2[j] - v0[j]
		}
		n := m.m_normals[i:]
		n[0] = e0[1]*e1[2] - e0[2]*e1[1]
		n[1] = e0[2]*e1[0] - e0[0]*e1[2]
		n[2] = e0[0]*e1[1] - e0[1]*e1[0]
		d := math.Sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2])
		if d > 0 {
			d = 1.0 / d
			n[0] *= d
			n[1] *= d
			n[2] *= d
		}
	}
}

func (m *rcMeshLoaderObj) addVertex(x, y, z float32) {
	m.m_verts = append(m.m_verts, x*m.m_scale, y*m.m_scale, z*m.m_scale)
	m.m_vertCount++
}

func (m *rcMeshLoaderObj) addTriangle(a, b, c int) {
	m.m_tris = append(m.m_tris, a, b, c)
	m.m_triCount++
}
