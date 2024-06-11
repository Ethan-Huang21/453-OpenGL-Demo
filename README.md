# OpenGL Demo
This is an academic project made from boiler-plate code provided during the 2021 Fall Semester for an Introductory Course to Computer Graphics (453) taken at the Universty of Calgary.
The project utilizes the OpenGL C++ Library and various third-party dependencies.
It features three key components:
 * An Interactive B-Spline & Bezier Curve Editor
 * A Camera-Movable Solid / Wireframe Model Viewer built from the Curve created in the Editor
 * A Hard-Coded Bezier Surface Viewer

# Execution
Run Cmake to get the necessary files -- Execute out/build/x64-Debug/453-skeleton.exe

# Controls
Mode Control (default: 0)
 * Right-Arrow: Increment Mode
 * Left-Arrow: Decrement Mode

Mode 0: Curve Editor
 * Alt: (default) [0] B-Spline | [1] Bezier
   * Up-Arrow: Increment Index
   * Down-Arrow: Decrement Index
 * Left-Click: Add Point | Hold+Drag Existing Point to Adjust
 * Right-Click: Remove Existing Point

Mode 1: Model Viewer
 * Alt: (default) [0] Curve Viewer | [1] Solid Model Viewer | [2] Wireframe Viewer
   * Up-Arrow: Increment Index
   * Down-Arrow: Decrement Index
 * Camera Control
   * W | A | S | D: Shift Camera Direction -- Forward | Backward | Left | Right
   * Mouse: Move to Look Around

Mode 2: Bezier Surface Viewer
 * Camera Control
   * W | A | S | D: Shift Camera Direction -- Forward | Backward | Left | Right
   * Mouse: Move to Look Around
   
