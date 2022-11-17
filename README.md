# UAV_stitching

This is An Ton's implementation to transfer all code to C++.

## How to use

The program is written with Visual Studio 2019 and OpenCV 4.5.5.

## Install OpenCV and others softwares to start

- Follow the [link](https://youtu.be/-GY2gT2umpk) to install OpenCV and OpenCV-contrib (please choose 4.5.5 version).

## Add environment variable's PATH

- After successfully install OpenCV and contrib, add *\<Your_folder>\opencv\build\install\x64\vc16\bin* to the PATH variable.

## Add directories into Visual Studio

- Create a new solution project on Visual Studio, and then add OpenCV's directories.
- Please change the **Solution configuration** from Debug -> Release, and the **Solution platform** into x64.
- Click to choose the solution project, then right-click **Project -> Properties**.
- In tab **Configuration Properties -> VC++ Directories**, at section **General -> Include Directories**, add *\<Your folder>\opencv\build\install\include*. At section **Library directory**, add *\<Your folder>\opencv\build\install\x64\vc16\lib*.
- In tab **Linker ->Input**, at section **Additional Dependencies**, include all *.lib* files in OpenCV directories *opencv\build\install\x64\vc16\lib*. It's usually that you can have only one `opencv_world4.x.x.lib` is enough.

## Code flow

In `main.cpp`

- At line 11, please specify the canvas directory path for the canvas image. The input images can be putted into folder `imgs` in the Github repo.
- At line 19 and 20, please specify the number of sub-images' center points and radiuses.
- At line 30, please specify the number of seconds and the fps of it. So the total frames after stitching together will be `frames = seconds x fps`.
- Also at line 30, please specify the changes in x and y coordinates after each frame to mimic multiple moving drones. *You can read more about the ImageSimulator class and its functions in `ImageSimulator.h` and `ImageSimulator.cpp`*.

Main flow: At the moment, after specify all above details, the code should generate an output video showing stitched images from 3 camera views and 30 seconds.

## Improvements

- The current implementation for 3 camera views (performing 2 stitches) takes around 6 to 8 seconds. However, to achieve real-time implementation, GPU can be thought of.
