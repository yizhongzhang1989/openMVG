#pragma once
#ifndef VISUALIZER_H_
#define VISUALIZER_H_

#include <vector>
#include <functional>
#include <mutex>
#include <thread>
#include <pangolin/pangolin.h>

namespace slam_visualization
{

using namespace pangolin;

template<class PointType, class PoseType, class AllocatorPoint = std::allocator<PointType>, class AllocatorPose = std::allocator<PoseType>>
class Visualizer
{
public:
    typedef PointType Point;
    typedef std::pair<Point, Point> Line;
    typedef PoseType Pose;

    typedef std::vector<PointType, AllocatorPoint> Points;
    typedef std::vector<Line, AllocatorPoint> Lines;
    typedef std::vector<PoseType, AllocatorPose> Poses;

    typedef std::function<OpenGlMatrix (const Pose&)> Pose2MatrixConverter;
    typedef std::function<Point (const Pose&)> PoseCenterGetter;

    struct VisualizerConfig
    {
        int width, height;
        double ViewpointX, ViewpointY, ViewpointZ, ViewpointF;
        double PointSize, LineWidth, KeyFrameSize, KeyFrameLineWidth, GraphLineWidth;
        double CameraSize, CameraLineWidth;
        Pose2MatrixConverter Pose2Matrix;
        PoseCenterGetter GetPoseCenter;

        VisualizerConfig()
        : width(1024), height(768), ViewpointX(0), ViewpointY(-0.7), ViewpointZ(-1.8), ViewpointF(500),
        PointSize(2.0), LineWidth(1.0),KeyFrameSize(0.05), KeyFrameLineWidth(1.0), GraphLineWidth(0.9),
        CameraSize(0.08), CameraLineWidth(3.0), Pose2Matrix(nullptr), GetPoseCenter(nullptr)
        {
            ;
        }
    };

    struct Color
    {
        GLfloat red, green, blue;
        Color()
        : red(0.0), green(0.0), blue(0.0)
        {
            ;
        }
        Color(GLfloat red, GLfloat green, GLfloat blue)
        : red(red), green(green), blue(blue)
        {
            ;
        }
    };

    struct ColorOptions
    {
        Color point_color;
        Color line_color;
        Color pose_color;
        Color graph_color;
        ColorOptions()
        : point_color(0.0,0.0,0.0), line_color(0.0,0.0,0.0), pose_color(0.0,0.0,1.0), graph_color(0.0,1.0,0.0)
        {
            ;
        }
    };
public:
    Visualizer(const std::string& windowName, const VisualizerConfig& config, bool multi_thread = false)
            : windowName(windowName), cfg(config), s_cam(nullptr), d_cam(nullptr),
            stop(false), mptDrawingThread(nullptr), mbMultiThread(multi_thread)
    {
        if(!mbMultiThread)
        {
            CreateWindow();
        }
    }
    ~Visualizer()
    {
        if(!mbMultiThread)
        {
            QuitAll();
        }
    }

    void CreateWindow()
    {
        CreateWindowAndBind(windowName, cfg.width, cfg.height);

        // 3D Mouse handler requires depth testing to be enabled
        glEnable(GL_DEPTH_TEST);

        // Issue specific OpenGl we might need
        glEnable (GL_BLEND);
        glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        // Define Camera Render Object (for view / scene browsing)
        s_cam = new OpenGlRenderState(
                ProjectionMatrix(cfg.width, cfg.height, cfg.ViewpointF, cfg.ViewpointF, cfg.width/2, cfg.height/2, 0.1, 1000),
                ModelViewLookAt(cfg.ViewpointX, cfg.ViewpointY, cfg.ViewpointZ, 0, 0, 0, 0, -1.0, 0)
        );

        // Add named OpenGL viewport to window and provide 3D Handler
        d_cam = &CreateDisplay()
                .SetBounds(0.0, 1.0, Attach::Pix(175), 1.0, -double(cfg.width)/double(cfg.height))
                .SetHandler(new Handler3D(*s_cam));
    }

    void DrawPoints(const Points& points, Color color = Color())
    {
        glPointSize(cfg.PointSize);
        glColor3f(color.red, color.green, color.blue);
        glBegin(GL_POINTS);

        for(const Point & p : points)
        {
            glVertex3f(p(0), p(1), p(2));
        }

        glEnd();
    }

    void DrawLines(const Lines& lines, Color color = Color())
    {
        glLineWidth(cfg.LineWidth);
        glColor3f(color.red, color.green, color.blue);
        glBegin(GL_LINES);

        for(const Line & l : lines)
        {
            const Point& sp = l.first;
            const Point& ep = l.second;
            glVertex3f(sp(0), sp(1), sp(2));
            glVertex3f(ep(0), ep(1), ep(2));
        }

        glEnd();
    }

    void DrawPoses(const Poses& poses, bool isKeyFrame = false, Color color = Color())
    {
        if(!cfg.Pose2Matrix)
            return;

        const float w = (isKeyFrame ? cfg.KeyFrameSize : cfg.CameraSize);
        const float h = w * 0.75;
        const float z = w * 0.6;

        for(const Pose & pose : poses)
        {
            glPushMatrix();
            OpenGlMatrix Twc = cfg.Pose2Matrix(pose);

#ifdef HAVE_GLES
            glMultMatrixf(Twc.m);
#else
            glMultMatrixd(Twc.m);
#endif

            if(isKeyFrame)
            {
                glLineWidth(cfg.KeyFrameLineWidth);
            }
            else
            {
                glLineWidth(cfg.CameraLineWidth);
            }
            glColor3f(color.red, color.green, color.blue);
            glBegin(GL_LINES);

            glVertex3f(0,0,0); glVertex3f(w,h,z);
            glVertex3f(0,0,0); glVertex3f(w,-h,z);
            glVertex3f(0,0,0); glVertex3f(-w,-h,z);
            glVertex3f(0,0,0); glVertex3f(-w,h,z);

            glVertex3f(w,h,z); glVertex3f(w,-h,z);
            glVertex3f(-w,h,z); glVertex3f(-w,-h,z);
            glVertex3f(-w,h,z); glVertex3f(w,h,z);
            glVertex3f(-w,-h,z); glVertex3f(w,-h,z);

            // draw axis
            const double axis_length = 0.05;
            glColor3f(1.0,0.0,0.0);
            glVertex3f(0,0,0); glVertex3f(axis_length,0,0);
            glColor3f(0.0,1.0,0.0);
            glVertex3f(0,0,0); glVertex3f(0,axis_length,0);
            glColor3f(0.0,0.0,1.0);
            glVertex3f(0,0,0); glVertex3f(0,0,axis_length);

            glEnd();

            glPopMatrix();
        }
    }

    void DrawGraph(const Poses& poses, Color color = Color())
    {
        if(!cfg.GetPoseCenter)
            return;

        glLineWidth(cfg.GraphLineWidth);
        glColor3f(color.red, color.green, color.blue);
        glBegin(GL_LINES);

        for(size_t i = 1; i < poses.size(); i++)
        {
            Point begin = cfg.GetPoseCenter(poses[i - 1]);
            Point end = cfg.GetPoseCenter(poses[i]);
            glVertex3f(begin(0), begin(1), begin(2));
            glVertex3f(end(0), end(1), end(2));
        }

        glEnd();
    }

    // Multi-thread functions
    void MultiThreadDrawing()
    {
        if(!mbMultiThread)
        {
            throw std::runtime_error("MultiThread drawing is not enabled.");
            return;
        }

        CreateWindow();
        stop = false;

        while(!stop)
        {
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//            if(!mPoses.empty())
//            {
//                s_cam->Follow(cfg.Pose2Matrix(mPoses.back()));
//            }
            d_cam->Activate(*s_cam);
            glClearColor(1.0f,1.0f,1.0f,1.0f);

            {
                std::unique_lock<std::mutex> lock(mMutexPoints);
                DrawPoints(mPoints, cfg_color.point_color);
            }
            {
                std::unique_lock<std::mutex> lock(mMutexLines);
                DrawLines(mLines, cfg_color.line_color);
            }
            {
                std::unique_lock<std::mutex> lock(mMutexPoses);
                DrawPoses(mPoses, false, cfg_color.pose_color);
                DrawGraph(mPoses, cfg_color.graph_color);
            }

            FinishFrame();
        }

        QuitAll();
    }
//    void StartMultiThread()
//    {
//        if(mptDrawingThread)
//            throw std::runtime_error("A thread has already been started.");
//
//        mptDrawingThread = new std::thread(&this->MultiThreadDrawing,this);
//    }
    void StopMultiThread()
    {
        stop = true;
//        mptDrawingThread->join();
//        mptDrawingThread = nullptr;
    }
    void MultiThreadUpdatePoints(const Points& points)
    {
        std::unique_lock<std::mutex> lock(mMutexPoints);
        mPoints = points;
    }
    void MultiThreadUpdateLines(const Lines& lines)
    {
        std::unique_lock<std::mutex> lock(mMutexLines);
        mLines = lines;
    }
    void MultiThreadUpdatePoses(const Poses& poses)
    {
        std::unique_lock<std::mutex> lock(mMutexPoses);
        mPoses = poses;
    }
    void UpdateColorOptions(const ColorOptions& options)
    {
        cfg_color = options;
    }
    pangolin::OpenGlRenderState* GetRenderState() const
    {
        return s_cam;
    }
    pangolin::View* GetView() const
    {
        return d_cam;
    }

private:
    VisualizerConfig cfg;
    pangolin::OpenGlRenderState* s_cam;
    pangolin::View* d_cam;
    std::string windowName;

    // for multi-thread
    Points mPoints;
    Lines mLines;
    Poses mPoses;
    ColorOptions cfg_color;
    bool stop;
    bool mbMultiThread;
    std::mutex mMutexPoints, mMutexLines, mMutexPoses;
    std::thread* mptDrawingThread;
};

}  // namespace slam_visualization

#endif  // VISUALIZER_H_