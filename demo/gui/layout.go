package gui

import (
	"fmt"
	"github.com/AllenDang/giu"
	"github.com/go-gl/gl/v4.1-core/gl"
	"gonavamesh/common"
	"math"
	"time"
)

// ui界面的管理和控制
type layout struct {
	ui            *Gui
	showLog       bool
	showMenu      bool
	showTools     bool
	showLevels    bool
	showSample    bool
	showTestCases bool

	// Window scroll positions.
	propScroll  int
	logScroll   int
	toolsScroll int

	 mouseOverMenu bool
	files map[string]string

	 geom *InputGeom
	 sample*Sample
	 test * TestCase

	lastTick time.Time

	cameraEulers []float64
	cameraPos  []float64
	camr float64
	meshesFolder string
	meshName string
}
func newLayout(ui *Gui) *layout {
	return &layout{
		meshesFolder : "Meshes",
		meshName : "Choose Mesh...",
		cameraEulers : []float64{45, -45},
		cameraPos :[]float64{0, 0, 0},
		camr : 1000,
		ui: ui,
		lastTick: time.Now(),
		showTools: true,
		showMenu: true}
}

func (l *layout) init() {

}

func (l *layout) Build() {
	var mouseScroll int
	 mousePos := []int{0, 0};
	 origMousePos := []int {0, 0}; // Used to compute mouse movement totals across frames.
	 origCameraEulers := []float64{0, 0}; // Used to compute rotational changes across frames.
	 moveFront := 0.0
	moveBack := 0.0
	moveLeft := 0.0
	moveRight := 0.0
	moveUp := 0.0
	moveDown := 0.0;
	  rayStart:=make([]float64,3)
	 rayEnd:=make([]float64,3)
	scrollZoom := 0
	timeAcc := 0.0
	done:=false
	 markerPosition := []float64{0, 0, 0};
	 markerPositionSet := false;
	var test *TestCase
	for done{
		now:=time.Now()

		// Handle input events.
		 mouseScroll := 0;
		var processHitTest bool
		var processHitTestShift bool

		dt:=now.Sub(l.lastTick).Seconds()
		// Hit test mesh.
		if (processHitTest && l.geom!=nil && l.sample!=nil) {
			var  hitTime float64
			 hit := l.geom.raycastMesh(rayStart, rayEnd, &hitTime);

			if (hit) {
				if (SDL_GetModState() & KMOD_CTRL) {
					// Marker
					markerPositionSet = true;
					markerPosition[0] = rayStart[0] + (rayEnd[0] - rayStart[0]) * hitTime;
					markerPosition[1] = rayStart[1] + (rayEnd[1] - rayStart[1]) * hitTime;
					markerPosition[2] = rayStart[2] + (rayEnd[2] - rayStart[2]) * hitTime;
				} else {
					 pos:=make([]float64,3)
					pos[0] = rayStart[0] + (rayEnd[0] - rayStart[0]) * hitTime;
					pos[1] = rayStart[1] + (rayEnd[1] - rayStart[1]) * hitTime;
					pos[2] = rayStart[2] + (rayEnd[2] - rayStart[2]) * hitTime;
					l.sample.handleClick(rayStart, pos, processHitTestShift);
				}
			} else {
				if (SDL_GetModState() & KMOD_CTRL) {
					// Marker
					markerPositionSet = false;
				}
			}
		}

		// Update sample simulation.
		 SIM_RATE := 20.0;
		 DELTA_TIME := 1.0 / SIM_RATE;
		timeAcc = common.Clamp(timeAcc + dt, -1.0, 1.0);
		 simIter := 0;
		for (timeAcc > DELTA_TIME){
			timeAcc -= DELTA_TIME;
			if (simIter < 5 && l.sample!=nil) {
				l.sample.handleUpdate(DELTA_TIME);
			}
			simIter++;
		}

		// Clamp the framerate so that we do not hog all the CPU.
		 MIN_FRAME_TIME := 1.0 / 40.0;
		if (dt < MIN_FRAME_TIME) {
			 ms := (int)((MIN_FRAME_TIME - dt) * 1000.0);
			if (ms > 10) {ms = 10;}
			if (ms >= 0) {gl.Delay(ms);}
		}

		// Set the viewport.
		glViewport(0, 0, l.ui.width, l.uiheight);
		GLint viewport[4];
		glGetIntegerv(GL_VIEWPORT, viewport);

		// Clear the screen
		glClearColor(0.3, 0.3, 0.32, 1.0);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		glDisable(GL_TEXTURE_2D);
		glEnable(GL_DEPTH_TEST);

		// Compute the projection matrix.
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(50.0, float64(l.ui.width)/float64(l.ui.height), 1.0, l.camr);
		GLdouble projectionMatrix[16];
		glGetDoublev(GL_PROJECTION_MATRIX, projectionMatrix);

		// Compute the modelview matrix.
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glRotatef(cameraEulers[0], 1, 0, 0);
		glRotatef(cameraEulers[1], 0, 1, 0);
		glTranslatef(-cameraPos[0], -l.cameraPos[1], -l.cameraPos[2]);
		GLdouble modelviewMatrix[16];
		glGetDoublev(GL_MODELVIEW_MATRIX, modelviewMatrix);

		// Get hit ray position and direction.
		var  x, y, z float64
		gluUnProject(mousePos[0], mousePos[1], 0.0, modelviewMatrix, projectionMatrix, viewport, &x, &y, &z);
		rayStart[0] = x;
		rayStart[1] = y;
		rayStart[2] = z;
		gluUnProject(mousePos[0], mousePos[1], 1.0, modelviewMatrix, projectionMatrix, viewport, &x, &y, &z);
		rayEnd[0] = x;
		rayEnd[1] = y;
		rayEnd[2] = z;

		// Handle keyboard movement.
		const Uint8* keystate = SDL_GetKeyboardState(NULL);
		moveFront	= rcClamp(moveFront	+ dt * 4 * ((keystate[SDL_SCANCODE_W] || keystate[SDL_SCANCODE_UP		]) ? 1 : -1), 0.0f, 1.0f);
		moveLeft	= rcClamp(moveLeft	+ dt * 4 * ((keystate[SDL_SCANCODE_A] || keystate[SDL_SCANCODE_LEFT		]) ? 1 : -1), 0.0f, 1.0f);
		moveBack	= rcClamp(moveBack	+ dt * 4 * ((keystate[SDL_SCANCODE_S] || keystate[SDL_SCANCODE_DOWN		]) ? 1 : -1), 0.0f, 1.0f);
		moveRight	= rcClamp(moveRight	+ dt * 4 * ((keystate[SDL_SCANCODE_D] || keystate[SDL_SCANCODE_RIGHT	]) ? 1 : -1), 0.0f, 1.0f);
		moveUp		= rcClamp(moveUp	+ dt * 4 * ((keystate[SDL_SCANCODE_Q] || keystate[SDL_SCANCODE_PAGEUP	]) ? 1 : -1), 0.0f, 1.0f);
		moveDown	= rcClamp(moveDown	+ dt * 4 * ((keystate[SDL_SCANCODE_E] || keystate[SDL_SCANCODE_PAGEDOWN	]) ? 1 : -1), 0.0f, 1.0f);

		 keybSpeed := 22.0;
		if (SDL_GetModState() & KMOD_SHIFT) {
			keybSpeed *= 4.0;
		}

		movex := (moveRight - moveLeft) * keybSpeed * dt;
		 movey := (moveBack - moveFront) * keybSpeed * dt + scrollZoom * 2.0;
		scrollZoom = 0;

		cameraPos[0] += movex *modelviewMatrix[0];
		cameraPos[1] += movex * modelviewMatrix[4];
		cameraPos[2] += movex * modelviewMatrix[8];

		cameraPos[0] += movey * modelviewMatrix[2];
		cameraPos[1] += movey * modelviewMatrix[6];
		cameraPos[2] += movey * modelviewMatrix[10];

		l.cameraPos[1] += (moveUp - moveDown) * keybSpeed * dt;

		glEnable(GL_FOG);

		if (l.sample!=nil){
			l.sample.handleRender();
		}

		if (l.test!=nil){
			l.test.handleRender();
		}


		glDisable(GL_FOG);

		// Render GUI
		glDisable(GL_DEPTH_TEST);
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluOrtho2D(0, width, 0, height);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();

		l.mouseOverMenu = false;
		point := giu.GetMousePos()
		l.ui.gs.beginFrame(point.X, point.Y, giu.IsMouseClicked(giu.MouseButtonLeft), mouseScroll)
		if (l.sample!=nil) {
			l.sample.handleRenderOverlay(projectionMatrix, modelviewMatrix, viewport);
		}
		if (test!=nil) {
			if (test.handleRenderOverlay(projectionMatrix, modelviewMatrix, viewport)){
			l.mouseOverMenu = true;
			}

		}
		l.ShowLevel()
		l.ShowTestCases()
		l.ShowMenu()
		l.ShowLog()
		l.ShowTools()
		// Marker
		if (markerPositionSet && gluProject(markerPosition[0], markerPosition[1], markerPosition[2],
			modelviewMatrix, projectionMatrix, viewport)){
		// Draw marker circle
		glLineWidth(5.0);
		glColor4ub(240,220,0,196);
		glBegin(GL_LINE_LOOP);
		 r := 25.0;
		for  i := 0; i < 20; i++{
		 a := float64(i) / 20.0 * math.Pi*2;
		 fx := x + math.Cos(a)*r;
		 fy :=y + math.Sin(a)*r;
		glVertex2f(fx,fy);
		}
		gl.End();
		gl.LineWidth(1.0);
		}
		l.ui.gs.endFrame()
		l.ui.render.Build()
		gl.Enable(gl.DEPTH_TEST);
		l.lastTick=now
	}


}

func (l *layout) ShowMenu() {
	 sampleName := "Choose Sample...";
	 testCasesFolder := "TestCases";
	// Help text.
	if !l.showMenu {
		return
	}
	msg := "W/S/A/D: Move  RMB: Rotate"
	l.ui.gs.imguiDrawText(280, l.ui.height-20, giu.AlignLeft, msg, imguiRGBA(255, 255, 255, 128))
	if l.ui.gs.imguiBeginScrollArea("Properties", l.ui.width-250-10, 10, 250, l.ui.height-20, &l.propScroll) {
		l.mouseOverMenu = true
	}

	if l.ui.gs.imguiCheck("Show Log", l.showLog, true) {
		l.showLog = !l.showLog
	}

	if l.ui.gs.imguiCheck("Show Tools", l.showTools, true) {
		l.showTools = !l.showTools
	}
	l.ui.gs.imguiSeparator()
	l.ui.gs.imguiLabel("Sample")
	if l.ui.gs.imguiButton(sampleName) {
		if l.showSample {
			l.showSample = false
		} else {
			l.showSample = true
			l.showLevels = false
			l.showTestCases = false
		}
	}

	l.ui.gs.imguiSeparator()
	l.ui.gs.imguiLabel("Input Mesh")
	if l.ui.gs.imguiButton(l.meshName) {
		if l.showLevels {
			l.showLevels = false
		} else {
			l.showSample = false
			l.showTestCases = false
			l.showLevels = true
			l.files= make(map[string]string)
			GetFileName(l.meshesFolder, ".obj",l.files)
			GetFileName(l.meshesFolder, ".gset",l.files)
		}
	}
	if l.geom !=nil{
		text := fmt.Sprintf("Verts: %.1fk  Tris: %.1fk",
			l.geom.getMesh().getVertCount()/1000.0,
			l.geom.getMesh().getTriCount()/1000.0)
		l.ui.gs.imguiValue(text)
	}
	l.ui.gs.imguiSeparator()

	if l.geom !=nil&& l.sample!=nil {
		l.ui.gs.imguiSeparatorLine()

		l.sample.handleSettings()

		if l.ui.gs.imguiButton("Build") {
			ctx.resetLog()
			if !l.ui.sample.handleBuild() {
				l.showLog = true
				l.logScroll = 0
			}
			ctx.dumpLog("Build log %s:", l.meshName)

			// Clear test.
			l.test = nil
		}

		l.ui.gs.imguiSeparator()
	}

	if l.ui.sample != nil {
		l.ui.gs.imguiSeparatorLine()
		l.ui.sample.handleDebugMode()
	}

	l.ui.gs.imguiEndScrollArea()
}

// Level selection dialog.
func (l *layout) ShowLevel() {
	 levelScroll := 0;
	if (l.ui.gs.imguiBeginScrollArea("Choose Level", l.ui.width - 10 - 250 - 10 - 200, l.ui.height - 10 - 450, 200, 450, &levelScroll)){
		l.mouseOverMenu = true;
	}

	exist:=false
	for k:=range l.files{
		if (l.ui.gs.imguiItem(k)) {
			l.meshName = k;
			exist=true
			break
	}
	}
	if (exist) {

		l.showLevels = false;
		l.geom = nil
		l.geom = newInputGeom()
		if (!l.geom.load( l.files[l.meshName])) {
			l.geom = nil
			// Destroy the sample if it already had geometry loaded, as we've just deleted it!
			if (l.sample!=nil && l.sample.getInputGeom()!=nil) {
				l.sample = nil
			}

			l.showLog = true;
			l.logScroll = 0;
			ctx.dumpLog("Geom load log %s:", l.meshName);
		}
		if (l.sample!=nil && l.geom!=nil) {
			l.sample.handleMeshChanged(l.geom);
		}

		if (l.geom!=nil || l.sample!=nil) {
			var bmin []float64
			var  bmax []float64
			if (l.geom!=nil) {
				bmin = l.geom.getNavMeshBoundsMin();
				bmax = l.geom.getNavMeshBoundsMax();
			}
			// Reset camera and fog to match the mesh bounds.
			if (bmin!=nil && bmax!=nil) {
				l.camr = math.Sqrt(common.Sqr(bmax[0]-bmin[0]) +
					common.Sqr(bmax[1]-bmin[1]) +
					common.Sqr(bmax[2]-bmin[2])) / 2;
				l.cameraPos[0] = (bmax[0] + bmin[0]) / 2 + l.camr;
				l.cameraPos[1] = (bmax[1] + bmin[1]) / 2 + l.camr;
				l.cameraPos[2] = (bmax[2] + bmin[2]) / 2 + l.camr;
				l.camr *= 3;
			}
			l.cameraEulers[0] = 45;
			l.cameraEulers[1] = -45;
			glFogf(GL_FOG_START, l.camr * 0.1);
			glFogf(GL_FOG_END, l.camr * 1.25);
		}
	}

	l.ui.gs.imguiEndScrollArea();
}

// Test cases
func (l *layout) ShowTestCases() {
	 testScroll := 0;
	if (l.ui.gs.imguiBeginScrollArea("Choose Test To Run", l.ui.width-10-250-10-200, l.ui.height-10-450, 200, 450, &testScroll)){
		l.mouseOverMenu = true;
	}

	path:=""
	for k,v:=range l.files{
		if l.ui.gs.imguiItem(k){
			path=v
		}
	}

	if (path!="") {
		l.test = newTestCase()
		if (l.test!=nil) {
			// Load the test.
			if (!l.test.load(path)) {
				l.test = nil
			}
			// Create sample
			var newSample *Sample
			for  i := 0; i < g_nsamples;i ++{
			if (g_samples[i].name == l.test.getSampleName()) {
			newSample = g_samples[i].create();
			if (newSample!=nil){
				sampleName = g_samples[i].name;
			}

			}
			}


			l.sample = newSample

			if (l.sample!=nil) {
				l.sample.setContext(&ctx);
				l.showSample = false;
			}

			// Load geom.
			l.meshName = l.test.getGeomFileName();


			path := l.meshesFolder + "/" + l.meshName;

			l.geom = newInputGeom()
			if (l.geom!=nil || !l.geom.load( path)) {
				l.geom = nil
				l.sample = nil
				l.showLog = true;
				l.logScroll = 0;
				ctx.dumpLog("Geom load log %s:", l.meshName);
			}
			if (l.sample!=nil && l.geom!=nil) {
				l.sample.handleMeshChanged(l.geom);
			}

			// This will ensure that tile & poly bits are updated in tiled sample.
			if (l.sample!=nil){
				l.sample.handleSettings();
			}


			ctx.resetLog();
			if (l.sample!=nil && !l.sample.handleBuild()) {
				ctx.dumpLog("Build log %s:", l.meshName);
			}

			if (l.geom!=nil || l.sample!=nil) {
				var  bmin []float64
				var  bmax []float64
				if (l.geom!=nil) {
					bmin =l.geom.getNavMeshBoundsMin();
					bmax = l.geom.getNavMeshBoundsMax();
				}
				// Reset camera and fog to match the mesh bounds.
				if (bmin!=nil && bmax!=nil) {
					l.camr = math.Sqrt(common.Sqr(bmax[0] - bmin[0]) +
						common.Sqr(bmax[1] - bmin[1]) +
						common.Sqr(bmax[2] - bmin[2])) / 2;
					l.cameraPos[0] = (bmax[0] + bmin[0]) / 2 + l.camr;
					l.cameraPos[1] = (bmax[1] + bmin[1]) / 2 + l.camr;
					l.cameraPos[2] = (bmax[2] + bmin[2]) / 2 + l.camr;
					l.camr *= 3;
				}
				l.cameraEulers[0] = 45;
				l.cameraEulers[1] = -45;
				glFogf(GL_FOG_START, l.camr * 0.2);
				glFogf(GL_FOG_END,l. camr * 1.25);
			}

			// Do the tests.
			if (l.sample!=nil){
				l.test.doTests(l.sample.getNavMesh(), l.sample.getNavMeshQuery());
			}

		}
	}

	l.ui.gs.imguiEndScrollArea();
}
func (l *layout) ShowSample(){
	 levelScroll := 0;
	if (l.ui.gs.imguiBeginScrollArea("Choose Sample", l.ui.width-10-250-10-200, l.ui.height-10-250, 200, 250, &levelScroll)){
		l.mouseOverMenu = true;
	}


	Sample* newSample = 0;
	for  i := 0; i < g_nsamples; i++{
	if (l.ui.gs.imguiItem(g_samples[i].name.c_str())) {
	newSample = g_samples[i].create();
	if (newSample)
	sampleName = g_samples[i].name;
	}
	}
	if (newSample) {
		l.sample = newSample;
		l.sample.setContext(&ctx);
		if (l.geom!=nil) {
			l.sample.handleMeshChanged(l.geom);
		}
		showSample = false;
	}

	if (l.geom!=nil || l.sample!=nil) {
		var bmin,bmax []float64
		if (l.geom!=nil) {
			bmin = l.geom.getNavMeshBoundsMin();
			bmax = l.geom.getNavMeshBoundsMax();
		}
		// Reset camera and fog to match the mesh bounds.
		if (bmin!=nil && bmax!=nil) {
			l.camr = math.Sqrt(common.Sqr(bmax[0]-bmin[0]) +
				common.Sqr(bmax[1]-bmin[1]) +
				common.Sqr(bmax[2]-bmin[2])) / 2;
			l.cameraPos[0] = (bmax[0] + bmin[0]) / 2 + l.camr;
			l.cameraPos[1] = (bmax[1] + bmin[1]) / 2 + l.camr;
			l.cameraPos[2] = (bmax[2] + bmin[2]) / 2 + l.camr;
			l.camr *= 3;
		}
		l.cameraEulers[0] = 45;
		l.cameraEulers[1] = -45;
		glFogf(GL_FOG_START, l.camr*0.1);
		glFogf(GL_FOG_END, l.camr*1.25);
	}

	l.ui.gs.imguiEndScrollArea();
}
// Log
func (l *layout) ShowLog() {
	if !(l.showLog && l.showMenu) {
		return
	}
	if l.ui.gs.imguiBeginScrollArea("Log", 250+20, 10, l.ui.width-300-250, 200, &l.logScroll) {
		l.mouseOverMenu = true
	}

	for i := 0; i < l.ui.logger.getLogCount(); i++ {
		l.ui.gs.imguiLabel(l.ui.logger.getLogText(i))
	}

	l.ui.gs.imguiEndScrollArea()
}

// Left column tools menu
func (l *layout) ShowTools() {
	// Left column tools menu
	if !(!l.showTestCases && l.showTools && l.showMenu) { // && geom && sample)
		return
	}
	if l.ui.gs.imguiBeginScrollArea("Tools", 10, 10, 250, l.ui.height-20, &l.toolsScroll) {
		l.mouseOverMenu = true
	}

	if l.ui.sample != nil {
		l.ui.sample.handleTools()
	}

	l.ui.gs.imguiEndScrollArea()
}
