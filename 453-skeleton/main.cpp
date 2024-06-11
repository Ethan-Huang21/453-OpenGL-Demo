#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <limits>
#include <functional>
#include <math.h>

#include "Geometry.h"
#include "GLDebug.h"
#include "Log.h"
#include "ShaderProgram.h"
#include "Shader.h"
#include "Texture.h"
#include "Window.h"

#include "glm/glm.hpp"
#include "glm/gtc/type_ptr.hpp"

#define PI 3.14159265

// Curve Type:
//	0 - B-Spline
//	1 - Bezier Curve
int curveType = 0;

// View Type:
//	0 - Curve Viewer
//	1 - Surface of Rev (Filled)
//	2 - Surface of Rev (WireFrame)
int viewType = 0;

// Model Type:
//	0 - Example 1 (Wireframe)
//	1 - Example 1 (Fill)
//	2 - Example 2 (Wireframe)
//	3 - Example 2 (Fill)
int modelType = 0;

// Mode:
//	0 - Curve Editor
//	1 - Viewing Model (3D)
//	2 - Tensor Model Surface
int mode = 0;

glm::mat4 viewMatrix(1.0f);
glm::mat4 projMatrix(1.0f);
glm::mat4 modelMatrix;
glm::mat4 mvpMatrix(1.0f);


// We gave this code in one of the tutorials, so leaving it here too
void updateGPUGeometry(GPU_Geometry &gpuGeom, CPU_Geometry const &cpuGeom) {
	gpuGeom.bind();
	gpuGeom.setVerts(cpuGeom.verts);
	gpuGeom.setCols(cpuGeom.cols);
}

CPU_Geometry BezierCurve(CPU_Geometry& cPoints, int degree) {
	CPU_Geometry bezierCurve = cPoints;
	CPU_Geometry storage;
	for (float u = 0; u < 1; u += 0.01) {
		for (int i = 1; i < degree; i++) {
			for (int j = 0; j < degree - i; j++) {
				bezierCurve.verts.at(j) = (1 - u) * bezierCurve.verts.at(j) + u * bezierCurve.verts.at(j+1);
			}
		}
		storage.verts.push_back(bezierCurve.verts.at(0));
	}
	return(storage);
}

glm::vec3 deCasteljau(CPU_Geometry& cPoints, int degree, float step) {
	CPU_Geometry bezierCurve = cPoints;
	glm::vec3 storage;
	for (int i = 1; i < degree; i++) {
		for (int j = 0; j < degree - i; j++) {
			bezierCurve.verts.at(j) = (1 - step) * bezierCurve.verts.at(j) + step * bezierCurve.verts.at(j + 1);
		}
	}
	storage = bezierCurve.verts.at(0);
	return(storage);
}

glm::vec3 sphereCoords(float phi, float theta, float radius) {
	float x = radius * sin(phi) * cos(theta);	// Sin phi, cos theta
	float y = radius * cos(phi);				// cos phi
	float z = radius * sin(phi) * sin(theta);	// Sin phi, sin theta
	return(glm::vec3(x, y, z));
}

CPU_Geometry sphere(float radius) {
	CPU_Geometry storage;
	float step = 0.10f;
	// i --> phi, j --> theta
	for (float i = 0; i < 3.1415; i += step) {
		for (float j = 0; j < 2 * (3.1415); j += step) {
			// B-Left
			storage.verts.push_back(sphereCoords(i, j, radius));
			// B-Right
			storage.verts.push_back(sphereCoords(i, j + step, radius));
			// T-Right
			storage.verts.push_back(sphereCoords(i + step, j + step, radius));
			// B-Left
			storage.verts.push_back(sphereCoords(i, j, radius));
			// T-Left
			storage.verts.push_back(sphereCoords(i + step, j, radius));
			// T-Right
			storage.verts.push_back(sphereCoords(i + step, j + step, radius));
		}
	}
	return(storage);
}

// Part 4
CPU_Geometry BezierSurface(CPU_Geometry& cPoints, int degree1, int degree2) {
	CPU_Geometry bezierControlPtsL;
	CPU_Geometry bezierControlPtsR;
	CPU_Geometry RvalBezierBL;
	CPU_Geometry RvalBezierBR;
	CPU_Geometry RvalBezierTL;
	CPU_Geometry RvalBezierTR;
	CPU_Geometry storage;
	glm::vec3 BotL, BotR, TopL, TopR;
	/*
	for (float u = 0; u < 1; u += 0.01) {
		for (float v = 0; v < 1; v += 0.01) {
			for (int i = 0; i < degree1; i++) {
				for (int j = 0; j < degree2; j++) {
					// takes in Pi,[0...d'] into some CPU_Geometry, in the verts.
					bezierControlPts.verts.push_back(cPoints.verts.at(i * degree1 + j));
				}
				// Take in one vertices from deCasteljau, at step v
				RvalBezier.verts.push_back(deCasteljau(bezierControlPts, degree1, v));
				// Clear bezierControlPts for the next row.
				bezierControlPts.verts.clear();
			}
			// Call deCasteljau again, at step u, using the RvalBezier vertices this time
			// Again, RvalBezier is of degree 3, as we go through 3 points for each row, 3 times, to get 1 point.
			// Meaning RvalBezier store 3 Rval points.
			storage.verts.push_back(deCasteljau(RvalBezier, degree1, u));
			// Clear for the next u-v.
			RvalBezier.verts.clear();
		}
	}
	*/
	
	for (float u = 0; u < 1; u += 0.1) {
		for (float v = 0; v < 1; v += 0.1) {
			// lopp needs to run 3x to get 3 points for ri.
			for (int i = 0; i < degree1; i++) {
				for (int j = 0; j < degree2; j++) {
					// takes in Pi,[0...d'] into some CPU_Geometry, in the verts.
					bezierControlPtsL.verts.push_back(cPoints.verts.at(i * degree1 + j));
				}

				// Calculate Rvals at four different points.
				RvalBezierBL.verts.push_back(deCasteljau(bezierControlPtsL, degree1, v));
				RvalBezierBR.verts.push_back(deCasteljau(bezierControlPtsL, degree1, v));
				RvalBezierTL.verts.push_back(deCasteljau(bezierControlPtsL, degree1, v+0.1));
				RvalBezierTR.verts.push_back(deCasteljau(bezierControlPtsL, degree1, v+0.1));

				// Clear bezierControlPts for the next row.
				bezierControlPtsL.verts.clear();
			}
			// Call deCasteljau again, at step u, using the RvalBezier vertices this time
			// Again, RvalBezier is of degree 3, as we go through 3 points for each row, 3 times, to get 1 point.
			// Meaning RvalBezier store 3 Rval points.
			// Calculate deCasteljau at the four R-Value Bezier points.
			BotL = deCasteljau(RvalBezierBL, degree1, u);
			BotR = deCasteljau(RvalBezierBR, degree1, u+0.1);
			TopL = deCasteljau(RvalBezierTL, degree1, u);
			TopR = deCasteljau(RvalBezierTR, degree1, u+0.1);

			storage.verts.push_back(BotL);
			storage.verts.push_back(TopL);
			storage.verts.push_back(TopR);
			storage.verts.push_back(BotL);
			storage.verts.push_back(BotR);
			storage.verts.push_back(TopR);
			// Clear for the next u-v.
			RvalBezierBL.verts.clear();
			RvalBezierBR.verts.clear();
			RvalBezierTL.verts.clear();
			RvalBezierTR.verts.clear();
		}
	}
	
	return(storage);
}

CPU_Geometry BSpline(CPU_Geometry& cPoints, int iteration) {
	CPU_Geometry bSpline;
	CPU_Geometry storage;
	float fval1 = 0.75f;
	float fval2 = 0.25f;

	int j = 0;
	if (iteration == 1) {
		// Return cPoints if our iteration is equal to 1 (for recursive purposes)
		return(cPoints);
	} else {
		// Store the first point in the new set of points.
		storage.verts.push_back(cPoints.verts.at(0));
		for (int i = 0; i < cPoints.verts.size() - 1; i++) {
			// Store the 1/4 and 3/4 new points, based on the size of our control points.
			storage.verts.push_back(fval1 * cPoints.verts.at(i) + fval2 * cPoints.verts.at(i+1));
			storage.verts.push_back(fval2 * cPoints.verts.at(i) + fval1 * cPoints.verts.at(i+1));
		}
		// Store the last point in the new set of points.
		storage.verts.push_back(cPoints.verts.at(cPoints.verts.size()-1));
		// Recursvely call BSpline to decrease the number of iterations (smoothing the curve out)
		bSpline = BSpline(storage, iteration - 1);
		return(bSpline);
	}
}

// UNUSED -- Tried to make one, didn't work out.
CPU_Geometry BSplineSurface(CPU_Geometry& cPoints, int iteration, int degree) {
	float fval = 0.75f;
	float fval2 = 0.25f;

	CPU_Geometry bSplineSection;
	CPU_Geometry bSplineCalculation;
	CPU_Geometry storage;
	int colAccess = 3;
	/*
	// Every single Column of Splines, ran through BSpline
	for (int i = 0; i < cPoints.verts.size(); i++) {
		bSplineSection.verts.push_back(cPoints.verts.at(i));
		if (bSplineSection.verts.size() == 3) {
			bSplineCalculation = BSpline(bSplineSection, iteration);
			for (int k = 0; k < bSplineCalculation.verts.size(); k++) {
				storage.verts.push_back(bSplineCalculation.verts.at(k));
			}
			bSplineSection.verts.clear();
		}
	}

	// Every single Row of Splines, ran through BSpline
	for (int i = 0; i < degree; i++) {
		for (int j = 0; j < degree; j++) {
			// 0 + 0 (pos 0,0), 0 + 3 (pos 1,0), 0 + 6 (pos, 2,0) etc --> Row Access
			bSplineSection.verts.push_back(cPoints.verts.at(i+(j*colAccess)));
			if (bSplineSection.verts.size() == 3) {
				bSplineCalculation = BSpline(bSplineSection, iteration);
				for (int k = 0; k < bSplineCalculation.verts.size(); k++) {
					storage.verts.push_back(bSplineCalculation.verts.at(k));
				}
				bSplineSection.verts.clear();
			}
		}
	}
	*/
	for (int i = 0; i < degree; i++) {
		for (int j = 0; j < degree; j++) {
			// takes in Pi,[0...d'] into some CPU_Geometry, in the verts.
			bSplineSection.verts.push_back(cPoints.verts.at(i * degree + j));
		}
		// Take in one vertices from deCasteljau, at step v
		bSplineCalculation = BSpline(bSplineSection, iteration);
		for (int k = 0; k < bSplineCalculation.verts.size(); k++) {
			storage.verts.push_back(bSplineCalculation.verts.at(k));
		}
		bSplineSection.verts.clear();
	}
	return(storage);
}

CPU_Geometry SurfaceRevolution(CPU_Geometry& curve, int iteration) {
	CPU_Geometry storage;

	for (float v = 0; v < 2 * PI; v += 0.1) {
		for (int i = 0; i < curve.verts.size()-1; i++) {
			// 6 to form a quad
			// bottom left
			storage.verts.push_back(glm::vec3(
				curve.verts.at(i).x * cos(v),
				curve.verts.at(i).y,
				curve.verts.at(i).x * sin(v)
			));
			// top left
			storage.verts.push_back(glm::vec3(
				curve.verts.at(i+1).x * cos(v),
				curve.verts.at(i+1).y,
				curve.verts.at(i+1).x * sin(v)
			));

			// top right
			storage.verts.push_back(glm::vec3(
				curve.verts.at(i + 1).x * cos(v + 0.1),
				curve.verts.at(i + 1).y,
				curve.verts.at(i + 1).x * sin(v + 0.1)
			));

			// bottom left
			storage.verts.push_back(glm::vec3(
				curve.verts.at(i).x * cos(v),
				curve.verts.at(i).y,
				curve.verts.at(i).x * sin(v)
			));
			// bottom right
			storage.verts.push_back(glm::vec3(
				curve.verts.at(i).x * cos(v + 0.1),
				curve.verts.at(i).y,
				curve.verts.at(i).x * sin(v + 0.1)
			));

			// top right
			storage.verts.push_back(glm::vec3(
				curve.verts.at(i + 1).x * cos(v + 0.1),
				curve.verts.at(i + 1).y,
				curve.verts.at(i + 1).x * sin(v + 0.1)
			));
		}
	}
	return(storage);
}

// EXAMPLE CALLBACKS
class Assignment3 : public CallbackInterface {

public:
	// Taken from CPSC453-MouseInput --> screenDim before shader as screenDim declared before shader below.
	Assignment3(int screenWidth, int screenHeight) :
		screenDim(screenWidth, screenHeight),
		leftMousePressed(false),
		rightMousePressed(false),
		leftMouseHeld(false),
		reset(false),
		mForward(false),
		mBackward(false),
		mLeft(false),
		mRight(false),
		mUp(false),
		mDown(false)
	{}

	virtual void keyCallback(int key, int scancode, int action, int mods) {
		if (key == GLFW_KEY_UP && action == GLFW_PRESS) {
			switch (mode) {
			case 0:
				curveType = 1;
				break;
			case 1:
				if (viewType < 2) {
					viewType++;
				}
				break;
			case 2:
				if (modelType < 3) {
					modelType++;
				}
				break;
			default:
				break;
			}
		}

		if (key == GLFW_KEY_DOWN && action == GLFW_PRESS) {
			switch (mode) {
			case 0:
				curveType = 0;
				break;
			case 1:
				if (viewType > 0) {
					viewType--;
				}
				break;
			case 2:
				if (modelType > 0) {
					modelType--;
				}
				break;
			default:
				break;
			}
		}

		if (key == GLFW_KEY_LEFT && action == GLFW_PRESS) {
			if (mode > 0) {
				mode--;
			}
		}

		if (key == GLFW_KEY_RIGHT && action == GLFW_PRESS) {
			if (mode < 2) {
				mode++;
			}
		}

		if (key == GLFW_KEY_R && action == GLFW_PRESS) {
			reset = true;
		}

		if (key == GLFW_KEY_W && action == GLFW_PRESS) {
			mForward = true;
		}

		if (key == GLFW_KEY_S && action == GLFW_PRESS) {
			mBackward = true;
		}

		if (key == GLFW_KEY_A && action == GLFW_PRESS) {
			mLeft = true;
		}

		if (key == GLFW_KEY_D && action == GLFW_PRESS) {
			mRight = true;
		}

		if (key == GLFW_KEY_Q && action == GLFW_PRESS) {
			mUp = true;
		}

		if (key == GLFW_KEY_E && action == GLFW_PRESS) {
			mDown = true;
		}
	}

	virtual void mouseButtonCallback(int button, int action, int mods) {
		// Create Point
		if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
			leftMouseHeld = true;
		}

		if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS) {
			rightMousePressed = true;
		}

		if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_RELEASE) {
			leftMousePressed = true;
			leftMouseHeld = false;
		}

	}

	virtual void cursorPosCallback(double xpos, double ypos) {
		// Saves the position to private variables --> Pixel coords, not openGL
		xScreenPos = xpos;
		yScreenPos = ypos;
	}

	virtual void scrollCallback(double xoffset, double yoffset) {
	}

	virtual void windowSizeCallback(int width, int height) {
		screenDim = glm::ivec2(width, height);
		// The CallbackInterface::windowSizeCallback will call glViewport for us
		CallbackInterface::windowSizeCallback(width,  height);
	}

	glm::vec3 mouseGL() {
		glm::vec2 startingVec(xScreenPos, yScreenPos);
		glm::vec2 shiftedVec = startingVec + glm::vec2(0.5f, 0.5f);
		// GLM interprets division as component-wise divison --> x1 / x2, y1 / y2
		glm::vec2 scaledToZeroOne = shiftedVec / glm::vec2(screenDim);

		glm::vec2 flippedY = glm::vec2(scaledToZeroOne.x, 1.0f - scaledToZeroOne.y);

		// 0-2 to -1 to 1.
		glm::vec2 final = flippedY * 2.0f - glm::vec2(1.0f, 1.0f);
		return(glm::vec3(final, 0.0f));
	}

	void refreshStatuses() {
		leftMousePressed = false;
		rightMousePressed = false;
		reset = false;
		mForward = false;
		mBackward = false;
		mLeft = false;
		mRight = false;
		mUp = false;
		mDown = false;
	}

	bool keyPressed(int button) {
		if (button == GLFW_KEY_R) return reset;
		if (button == GLFW_KEY_W) return mForward;
		if (button == GLFW_KEY_S) return mBackward;
		if (button == GLFW_KEY_A) return mLeft;
		if (button == GLFW_KEY_D) return mRight;
		if (button == GLFW_KEY_Q) return mUp;
		if (button == GLFW_KEY_E) return mDown;
		return false;
	}

	bool mouseButtonDown(int button) {
		if (button == GLFW_MOUSE_BUTTON_LEFT) return leftMousePressed;
		else if (button == GLFW_MOUSE_BUTTON_RIGHT) return rightMousePressed;

		return false;
	}

	bool mouseButtonHeld(int button) {
		if (button == GLFW_MOUSE_BUTTON_LEFT) return leftMouseHeld;
		else return false;
	}

	bool pointColliding(glm::vec3 mousePos, CPU_Geometry& cPoints) {
		if (cPoints.verts.size() > 0) {
			for (int i = 0; i < cPoints.verts.size(); i++) {
				if (mousePos.x <= cPoints.verts.at(i).x + 0.08f
					&& mousePos.x >= cPoints.verts.at(i).x - 0.08f
					&& mousePos.y <= cPoints.verts.at(i).y + 0.08f
					&& mousePos.y >= cPoints.verts.at(i).y - 0.08f
					) {
					return(true);
				}
			}
			return(false);
		}
		else {
			return(false);
		}
	}

	int getCollidingPoint(glm::vec3 mousePos, CPU_Geometry& cPoints) {
		if (cPoints.verts.size() > 0) {
			for (int i = 0; i < cPoints.verts.size(); i++) {
				if (mousePos.x <= cPoints.verts.at(i).x + 0.08f
					&& mousePos.x >= cPoints.verts.at(i).x - 0.08f
					&& mousePos.y <= cPoints.verts.at(i).y + 0.08f
					&& mousePos.y >= cPoints.verts.at(i).y - 0.08f
					) {
					return(i);
				}
			}
			return(-1);
		}
		else {
			return(-1);
		}
	}

private:
	glm::ivec2 screenDim;

	double xScreenPos = 0.0f;
	double yScreenPos = 0.0f;

	// Can be made into some dictionary instead of booleans.
	bool leftMousePressed;
	bool rightMousePressed;
	bool leftMouseHeld;
	bool pointCollision;
	bool reset;
	bool mForward;
	bool mBackward;
	bool mLeft;
	bool mRight;
	bool mUp;
	bool mDown;
};

int main() {
	Log::debug("Starting main");
	int index = -1;
	// WINDOW
	glfwInit();
	Window window(800, 800, "CPSC 453"); // can set callbacks at construction if desired

	GLDebug::enable();

	// CALLBACKS
	auto a3 = std::make_shared<Assignment3>(800, 800);
	window.setCallbacks(a3);

	ShaderProgram shader("shaders/test.vert", "shaders/test.frag");

	// The current CPU_Geometry and GPU_Geometry classes are defined in
	// Geometry.h/Geometry.cpp They will work for this assignment, but for some of
	// the bonuses you may have to modify them.
	CPU_Geometry square;
	GPU_Geometry pointsGPUGeom;
	GPU_Geometry linesGPUGeom;
	glPointSize(10.0f);

	// CPU Geometry for the curve
	// GPU Geometry for the curve
	CPU_Geometry curveCPUGeom;
	GPU_Geometry curveGPUGeom;

	// CPU Geometry for the Surface of Revolution
	// GPU Geomery for the Surface of Revolution
	CPU_Geometry surfaceRevGeom;
	GPU_Geometry surfaceRevGPU;

	// CPU Geometry for the Tensor surfaces
	// GPU Geometry for the Tensor surfaces
	CPU_Geometry tensorSurfaceCPU;
	GPU_Geometry tensorSurfaceGPU;

	// Parametric Surface Vector of Vectors
	std::vector<glm::vec3> paramSurface;

	// Note this call only work on some systems, unfortunately.
	// In order for them to work, you have to comment out line 60
	// If you're on a mac, you can't comment out line 60, so you
	// these will have no effect. :(
	// glLineWidth(5.0f);

	// Camera Vectors
	bool firstMouse = true;
	glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 3.0f);
	glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
	glm::vec3 cameraUp = glm::vec3(0.0f, 1.0f, 0.0f);
	glm::vec3 direction;
	float lastX = 0;
	float lastY = 0;
	float yaw = -90.0f;
	float pitch = 0.0f;

	//CPU_Geometry SPH = sphere(1.0f);
	//GPU_Geometry SPHGPU;
	//updateGPUGeometry(SPHGPU, SPH);

	// RENDER LOOP
	while (!window.shouldClose()) {
		a3->refreshStatuses();
		glfwPollEvents();

		// Change points
		if (mode == 0) {
			// Clear the board
			if (a3->keyPressed(GLFW_KEY_R)) {
				square.verts.clear();
				square.cols.clear();
				curveCPUGeom.verts.clear();
				curveCPUGeom.cols.clear();
			}

			// Delete points by checking if hte right mouse button is clicked while colliding with existing points.
			if (a3->mouseButtonDown(GLFW_MOUSE_BUTTON_RIGHT) && a3->pointColliding(a3->mouseGL(), square)) {
				index = a3->getCollidingPoint(glm::vec3(a3->mouseGL()), square);
				square.verts.erase(square.verts.begin() + index);

				// Remove curve if there are less than 3 points.
				if (square.verts.size() < 3) {
					curveCPUGeom.verts.clear();
					curveCPUGeom.cols.clear();
				}
			}

			// Modify Points by checking of the left mouse button is held while colliding with existing points
			if (a3->mouseButtonHeld(GLFW_MOUSE_BUTTON_LEFT) && a3->pointColliding(a3->mouseGL(), square)) {
				index = a3->getCollidingPoint(glm::vec3(a3->mouseGL()), square);
				square.verts.at(index) = a3->mouseGL();
			}

			// Add points by checking if the left mouse button was clicked, while not colliding with any existing points
			if (a3->mouseButtonDown(GLFW_MOUSE_BUTTON_LEFT) && !(a3->pointColliding(a3->mouseGL(), square))) {
				square.verts.push_back(a3->mouseGL());
			}

			// Update the Curve
			if (square.verts.size() > 2) {
				switch (curveType) {
				case 0:
					curveCPUGeom = BSpline(square, 5);
					break;
				case 1:
					curveCPUGeom = BezierCurve(square, square.verts.size());
					break;
				default:
					curveCPUGeom = BSpline(square, 5);
					break;
				}
				curveCPUGeom.cols.clear();
				curveCPUGeom.cols.resize(curveCPUGeom.verts.size(), glm::vec3{ 0.0, 0.0, 1.0 });
				updateGPUGeometry(curveGPUGeom, curveCPUGeom);
			}

			projMatrix = glm::mat4(1.0f);

			// Reset Viewmatrix back to original view
			viewMatrix = glm::lookAt(
				glm::vec3(0.0f, 0.0f, 1.0f),
				glm::vec3(0.0f, 0.0f, 0.0f),
				glm::vec3(0.0f, 1.0f, 0.0f)
			);

			mvpMatrix = projMatrix * viewMatrix;
		}

		// Surface of Revolution & 3D Curve Viewer
		if (mode == 1) {
			// If the viewType is > 0 (meaning we're viewing the wireframe)
			// Switch the curveType back to B-Spline.
			if (curveType != 0 && viewType > 0) {
				curveType = 0;
				curveCPUGeom = BSpline(square, 5);
				curveCPUGeom.cols.clear();
				curveCPUGeom.cols.resize(curveCPUGeom.verts.size(), glm::vec3{ 0.0, 0.0, 1.0 });
				updateGPUGeometry(curveGPUGeom, curveCPUGeom);
			}

			// Surface of Revolution
			switch (viewType) {
			case 1:
				surfaceRevGeom = SurfaceRevolution(curveCPUGeom, 1);
				surfaceRevGeom.cols.clear();
				surfaceRevGeom.cols.resize(surfaceRevGeom.verts.size(), glm::vec3{ 0, 0, 0 });
				updateGPUGeometry(surfaceRevGPU, surfaceRevGeom);
				break;
			case 2:
				surfaceRevGeom = SurfaceRevolution(curveCPUGeom, 1);
				surfaceRevGeom.cols.clear();
				surfaceRevGeom.cols.resize(surfaceRevGeom.verts.size(), glm::vec3{ 0, 0, 0 });
				updateGPUGeometry(surfaceRevGPU, surfaceRevGeom);
				break;
			default:
				break;
			}

			//float camX = sin(glfwGetTime());
			//float camZ = cos(glfwGetTime());
			//glm::vec3(camX, 0.0f, camZ),
		}

		// Tensor Product Surfaces
		if (mode == 2) {
			switch (modelType) {
			case 0:
				square.verts.clear();
				square.cols.clear();
				square.verts.push_back(glm::vec3(0.0f, -0.5f, 0.0f));
				square.verts.push_back(glm::vec3(0.5f, -0.5f, 0.0f));
				square.verts.push_back(glm::vec3(1.0f, -0.5f, 0.0f));
				square.verts.push_back(glm::vec3(0.0f, 0.0f, 0.5f));
				square.verts.push_back(glm::vec3(0.5f, -0.25f, 0.5f));
				square.verts.push_back(glm::vec3(1.0f, -0.75f, 0.5f));
				square.verts.push_back(glm::vec3(0.0f, -0.5f, 1.0f));
				square.verts.push_back(glm::vec3(0.5f, -0.5f, 1.0f));
				square.verts.push_back(glm::vec3(1.0f, -0.5f, 1.0f));
				break;
			case 1:
				square.verts.clear();
				square.cols.clear();
				square.verts.push_back(glm::vec3(0.0f, -0.5f, 0.0f));
				square.verts.push_back(glm::vec3(0.5f, -0.5f, 0.0f));
				square.verts.push_back(glm::vec3(1.0f, -0.5f, 0.0f));
				square.verts.push_back(glm::vec3(0.0f, 0.0f, 0.5f));
				square.verts.push_back(glm::vec3(0.5f, -0.25f, 0.5f));
				square.verts.push_back(glm::vec3(1.0f, -0.75f, 0.5f));
				square.verts.push_back(glm::vec3(0.0f, -0.5f, 1.0f));
				square.verts.push_back(glm::vec3(0.5f, -0.5f, 1.0f));
				square.verts.push_back(glm::vec3(1.0f, -0.5f, 1.0f));
				break;
			case 2:
				square.verts.clear();
				square.cols.clear();
				square.verts.push_back(glm::vec3(0.0f, -0.5f, 0.0f));
				square.verts.push_back(glm::vec3(0.5f, 0.0f, 0.0f));
				square.verts.push_back(glm::vec3(1.0f, -0.5f, 0.0f));
				square.verts.push_back(glm::vec3(0.0f, 0.0f, 0.5f));
				square.verts.push_back(glm::vec3(0.5f, -0.5f, 0.5f));
				square.verts.push_back(glm::vec3(1.0f, 0.5f, 0.5f));
				square.verts.push_back(glm::vec3(0.0f, -0.5f, 1.0f));
				square.verts.push_back(glm::vec3(0.5f, -0.75f, 1.0f));
				square.verts.push_back(glm::vec3(1.0f, -0.5f, 1.0f));
				break;
			case 3:
				square.verts.clear();
				square.cols.clear();
				square.verts.push_back(glm::vec3(0.0f, -0.5f, 0.0f));
				square.verts.push_back(glm::vec3(0.5f, 0.0f, 0.0f));
				square.verts.push_back(glm::vec3(1.0f, -0.5f, 0.0f));
				square.verts.push_back(glm::vec3(0.0f, 0.0f, 0.5f));
				square.verts.push_back(glm::vec3(0.5f, -0.5f, 0.5f));
				square.verts.push_back(glm::vec3(1.0f, 0.5f, 0.5f));
				square.verts.push_back(glm::vec3(0.0f, -0.5f, 1.0f));
				square.verts.push_back(glm::vec3(0.5f, -0.75f, 1.0f));
				square.verts.push_back(glm::vec3(1.0f, -0.5f, 1.0f));
				break;
			default:
				break;
			}

			// Call Bezier Surface to get the CPU Geometry
			tensorSurfaceCPU = BezierSurface(square, 3, 3);
			// Set color to Black and Update GPU Geometry
			tensorSurfaceCPU.cols.clear();
			tensorSurfaceCPU.cols.resize(tensorSurfaceCPU.verts.size(), glm::vec3{ 0, 0, 0 });
			updateGPUGeometry(tensorSurfaceGPU, tensorSurfaceCPU);

		}

		// Camera Controls, for modes 1 and 2
		if (mode == 1 || mode == 2) {
			// REFERENCE:
			// https://learnopengl.com/Getting-started/Camera?fbclid=IwAR1DNMFvCq0ieQgH__qG1i7-BiiPvDyIPBrVF5xW4AwQhC7aoEgamFB2Y9Q
			// https://learnopengl.com/Getting-started/Coordinate-Systems?fbclid=IwAR1NVxwAFDNCtIpEeTXoRSyXxbtHxkOc9Pd460GTvSHWG-2UhjPW2eWh6Ns

			// Key Controls
			const float cameraSpeed = 0.05f;
			if (a3->keyPressed(GLFW_KEY_W)) {
				// Move Camera Position Forward
				cameraPos += cameraSpeed * cameraFront;
			}
			if (a3->keyPressed(GLFW_KEY_S)) {
				// Move Camera Position Backward
				cameraPos -= cameraSpeed * cameraFront;
			}
			if (a3->keyPressed(GLFW_KEY_A)) {
				// Move Camera Position Leftward
				cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
			}
			if (a3->keyPressed(GLFW_KEY_D)) {
				// Move Camera Position Rightward
				cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
			}
			if (a3->keyPressed(GLFW_KEY_Q)) {
				// Move Camera Position Upward
				cameraPos += cameraUp * cameraSpeed;
			}
			if (a3->keyPressed(GLFW_KEY_E)) {
				// Move Camera Position Downward
				cameraPos -= cameraUp * cameraSpeed;
			}

			// Mouse Control (TAKEN FROM REFERENCE ABOVE)
			if (a3->mouseButtonHeld(GLFW_MOUSE_BUTTON_LEFT)) {
				glm::vec3 mousePosition = a3->mouseGL();

				if (firstMouse) {
					lastX = mousePosition.x;
					lastY = mousePosition.y;
					firstMouse = false;
				}

				float xoffset = mousePosition.x - lastX;
				float yoffset = mousePosition.y - lastY;

				const float sensitivity = 0.2f;
				xoffset *= sensitivity;
				yoffset *= sensitivity;

				yaw += xoffset;
				pitch += yoffset;

				if (pitch > 89.0f) {
					pitch = 89.0f;
				}
				if (pitch < -89.0f) {
					pitch = -89.0f;
				}

				direction.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
				direction.y = sin(glm::radians(pitch));
				direction.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
				cameraFront = glm::normalize(direction);
			}

			// Camera Position, Camera Target, up Vector
			viewMatrix = glm::lookAt(
				cameraPos,
				cameraPos + cameraFront,
				cameraUp
			);

			// Multiply projMatrix and viewMatrix to create the matrix to pass to shaders
			mvpMatrix = projMatrix * viewMatrix;
		}

		// Send global variable 'viewMatrix' to the 'viewMatrix' mat4 in shaders.
		GLint myLoc = glGetUniformLocation(shader.getProgram(), "viewMatrix");
		glUniformMatrix4fv(myLoc, 1, GL_FALSE, &mvpMatrix[0][0]);

		// Perspective
		projMatrix = glm::perspective(glm::radians(45.0f), (float)800 / (float)800, 0.1f, 100.0f);

		if (square.verts.size() != 0) {
			// Set to Red for Points
			square.cols.clear();
			square.cols.resize(square.verts.size(), glm::vec3{ 1.0, 0.0, 0.0 });
			updateGPUGeometry(pointsGPUGeom, square);

			// Reset to green for Lines
			square.cols.clear();
			square.cols.resize(square.verts.size(), glm::vec3{ 0.0, 1.0, 0.0 });
			updateGPUGeometry(linesGPUGeom, square);
		}

		glEnable(GL_LINE_SMOOTH);
		glEnable(GL_FRAMEBUFFER_SRGB);
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		// switch polygonMode depending on viewType.
		// 1 means fill
		// 2 means wireframe
		switch (viewType) {
		case 1:
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
			break;
		case 2:
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			break;
		default:
			break;
		}
		if (mode == 2) {
			switch (modelType) {
			case 0:
				glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
				break;
			case 1:
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				break;
			case 2:
				glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
				break;
			case 3:
				glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
				break;
			default:
				break;
			}
		}

		shader.use();

		// Draw, Depending on which 'mode' we're on.
		if (mode == 0 && square.verts.size() != 0) {
			linesGPUGeom.bind();
			glDrawArrays(GL_LINE_STRIP, 0, GLsizei(square.verts.size()));

			curveGPUGeom.bind();
			glDrawArrays(GL_LINE_STRIP, 0, GLsizei(curveCPUGeom.verts.size()));

			pointsGPUGeom.bind();
			glDrawArrays(GL_POINTS, 0, GLsizei(square.verts.size()));
		}
		else if (mode == 1) {
			switch (viewType) {
			case 0:
				linesGPUGeom.bind();
				glDrawArrays(GL_LINE_STRIP, 0, GLsizei(square.verts.size()));

				curveGPUGeom.bind();
				glDrawArrays(GL_LINE_STRIP, 0, GLsizei(curveCPUGeom.verts.size()));

				pointsGPUGeom.bind();
				glDrawArrays(GL_POINTS, 0, GLsizei(square.verts.size()));
				break;
			case 1:
				// Fill Mode for Surface of Revolution
				surfaceRevGPU.bind();
				glDrawArrays(GL_TRIANGLES, 0, GLsizei(surfaceRevGeom.verts.size()));
				break;
			case 2:
				// Wireframe Mode for Surface of Revolution
				surfaceRevGPU.bind();
				glDrawArrays(GL_TRIANGLES, 0, GLsizei(surfaceRevGeom.verts.size()));
				break;
			default:
				break;
			}
		}
		else if (mode == 2) {
			pointsGPUGeom.bind();
			glDrawArrays(GL_POINTS, 0, GLsizei(square.verts.size()));

			tensorSurfaceGPU.bind();
			glDrawArrays(GL_TRIANGLES, 0, GLsizei(tensorSurfaceCPU.verts.size()));
		}

		glDisable(GL_FRAMEBUFFER_SRGB); // disable sRGB for things like imgui

		window.swapBuffers();
	}

	glfwTerminate();
	return 0;
}
