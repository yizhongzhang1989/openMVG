#include <iostream>
#include <thread>

#include "types.h"
#include "Visualizer.h"
#include "TrajectorySampler.h"

#ifdef USE_PANGOLIN
#include "Visualizer.h"
typedef slam_visualization::Visualizer<Eigen::Vector3d, generator::InversePose, Eigen::aligned_allocator<Eigen::Vector3d>, Eigen::aligned_allocator<generator::InversePose>> GeneratorVisualizer;
#endif

using namespace generator;
using namespace std;

void DrawGrid(STLVector<Eigen::Vector3d>& up_vertices, STLVector<Eigen::Vector3d>& down_vertices);

int main(int argc, char** argv)
{
    if(argc != 2)
    {
        cerr<<"usage: "<<argv[0]<<" trajectory.obj"<<endl;
        return EXIT_FAILURE;
    }

    string objFile(argv[1]);
    STLVector<Eigen::Vector3d> vertices;
    STLVector<InversePose> inv_poses = TrajectorySampler::SampleTrajectory(objFile,&vertices);

    int N = vertices.size() / 2;
    STLVector<Eigen::Vector3d> up_vertices, down_vertices;
    for(int i = 0; i < N; i++)
    {
        down_vertices.push_back(vertices[2 * i]);
        up_vertices.push_back(vertices[2 * i + 1]);
    }

#ifdef USE_PANGOLIN
    auto GetPoseCenter = [](const InversePose& pose)->Eigen::Vector3d
    {
        return pose.p;
    };
    auto Pose2Matrix = [](const InversePose& pose)->pangolin::OpenGlMatrix
    {
        Eigen::Matrix4d Twc;
        Twc << pose.q.toRotationMatrix(), pose.p, Eigen::Vector3d::Zero().transpose(), 1;
        return {Twc};
    };
    GeneratorVisualizer::VisualizerConfig cfg;
    cfg.GetPoseCenter = GetPoseCenter;
    cfg.Pose2Matrix = Pose2Matrix;
    GeneratorVisualizer visualizer("Visualization",cfg);

    pangolin::OpenGlRenderState* s_cam = visualizer.GetRenderState();
    pangolin::View* d_cam = visualizer.GetView();
    GeneratorVisualizer::ColorOptions color_opt;

    while(!pangolin::ShouldQuit())
    {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        d_cam->Activate(*s_cam);
        glClearColor(1.0f,1.0f,1.0f,1.0f);

        visualizer.DrawPoses(inv_poses,false,color_opt.pose_color);
        visualizer.DrawGraph(inv_poses,color_opt.graph_color);
        DrawGrid(up_vertices,down_vertices);

        pangolin::FinishFrame();
    }

#endif

    return 0;
}

void DrawGrid(STLVector<Eigen::Vector3d>& up_vertices, STLVector<Eigen::Vector3d>& down_vertices)
{
    assert(up_vertices.size() == down_vertices.size());
    int N = up_vertices.size();

    glLineWidth(2);
    glBegin(GL_LINES);

    for(int i = 0; i < N - 1; i++)
    {
        Eigen::Vector3d& O = down_vertices[i];
        Eigen::Vector3d& A = down_vertices[i+1];
        Eigen::Vector3d& B = up_vertices[i+1];
        Eigen::Vector3d& C = up_vertices[i];

        glColor3f(1.0,0,1.0);
        glVertex3f(O[0],O[1],O[2]); glVertex3f(A[0],A[1],A[2]);
        glVertex3f(C[0],C[1],C[2]); glVertex3f(B[0],B[1],B[2]);

        glColor3f(0,1.0,1.0);
        glVertex3f(O[0],O[1],O[2]); glVertex3f(C[0],C[1],C[2]);
        glVertex3f(A[0],A[1],A[2]); glVertex3f(B[0],B[1],B[2]);
    }

    glEnd();
}