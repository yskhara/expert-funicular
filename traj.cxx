//////////////////////////////////////////////////////////////////////////
//////////////////            lts.cxx        /////////////////////////////
//////////////////////////////////////////////////////////////////////////
////////////////           PSOPT  Example             ////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////// Title: Linear tangent steering problem           ////////////////
//////// Last modified:         16 February 2009          ////////////////
//////// Reference:             Betts (2001)              ////////////////
//////// (See PSOPT handbook for full reference)          ////////////////
//////////////////////////////////////////////////////////////////////////
////////     Copyright (c) Victor M. Becerra, 2009        ////////////////
//////////////////////////////////////////////////////////////////////////
//////// This is part of the PSOPT software library, which ///////////////
//////// is distributed under the terms of the GNU Lesser ////////////////
//////// General Public License (LGPL)                    ////////////////
//////////////////////////////////////////////////////////////////////////

#include "psopt.h"
#include <GLFW/glfw3.h>

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the end point (Mayer) cost function //////////
//////////////////////////////////////////////////////////////////////////

adouble endpoint_cost(adouble* initial_states, adouble* final_states, adouble* parameters, adouble& t0, adouble& tf,
        adouble* xad, int iphase, Workspace* workspace)
{
    return tf;
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the integrand (Lagrange) cost function  //////
//////////////////////////////////////////////////////////////////////////

adouble integrand_cost(adouble* states, adouble* controls, adouble* parameters, adouble& time, adouble* xad, int iphase,
        Workspace* workspace)
{
    double w = 1.0;

    return 0;//w * (pow(controls[CINDEX(1)], 2) + pow(controls[CINDEX(2)], 2));
}

//////////////////////////////////////////////////////////////////////////
///////////////////  Define the DAE's ////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

void dae(adouble* derivatives, adouble* path, adouble* states, adouble* controls, adouble* parameters, adouble& time,
        adouble* xad, int iphase, Workspace* workspace)
{

    adouble x1 = states[CINDEX(1)];     // world x
    adouble x2 = states[CINDEX(2)];     // world y
    adouble x3 = states[CINDEX(3)];     // theta
    adouble x4 = states[CINDEX(4)];     // d/dt phi_1
    adouble x5 = states[CINDEX(5)];     // d/dt phi_2
    //adouble x6 = states[CINDEX(6)];

    adouble u1 = controls[CINDEX(1)];   // u_1, torque on wheel 1
    adouble u2 = controls[CINDEX(2)];   // u_2, torque on wheel 2

    double m = 10;                      // mass of the robot (incl. wheels)
    double L = 0.1;
    double r = 0.05;
    double i_z = 0.1 * r * r;
    double i_y = 0.5 * r * r / 2;       // inertia of the wheels

    derivatives[CINDEX(1)] = r * ((x4 * cos(x3)) - (x5 * sin(x3))) / 2;
    derivatives[CINDEX(2)] = r * ((x4 * sin(x3)) + (x5 * cos(x3))) / 2;
    derivatives[CINDEX(3)] = L * (x4 - x5) / (2 * r);

    adouble trans_force = (u1 + u2) / (m * r * r + 2 * i_y);
    adouble rot_force = L * L * (u1 - u2) / (i_z * r * r + 2 * L * L * i_y);

    derivatives[CINDEX(4)] = trans_force + rot_force;
    derivatives[CINDEX(5)] = trans_force - rot_force;
    //derivatives[CINDEX(6)] = x3;

    path[CINDEX(1)] = L * (x4 - x5) / (2 * r); // d/dt theta
    path[CINDEX(2)] = pow(x1 - 1, 2.0) + pow(x2 - 1, 2.0);
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the events function ////////////////////////////
////////////////////////////////////////////////////////////////////////////

void events(adouble* e, adouble* initial_states, adouble* final_states, adouble* parameters, adouble& t0, adouble& tf,
        adouble* xad, int iphase, Workspace* workspace)
{
    adouble x10 = initial_states[CINDEX(1)];
    adouble x20 = initial_states[CINDEX(2)];
    adouble x30 = initial_states[CINDEX(3)];
    adouble x40 = initial_states[CINDEX(4)];
    adouble x50 = initial_states[CINDEX(5)];
    //adouble x60 = initial_states[CINDEX(6)];
    adouble x1f = final_states[CINDEX(1)];
    adouble x2f = final_states[CINDEX(2)];
    adouble x3f = final_states[CINDEX(3)];
    adouble x4f = final_states[CINDEX(4)];
    adouble x5f = final_states[CINDEX(5)];
    //adouble x6f = final_states[CINDEX(6)];

    e[CINDEX(1)] = x10;
    e[CINDEX(2)] = x20;
    e[CINDEX(3)] = x30;
    e[CINDEX(4)] = x40;
    e[CINDEX(5)] = x50;
    e[CINDEX(6)] = x1f;
    e[CINDEX(7)] = x2f;
    e[CINDEX(8)] = x3f;
    e[CINDEX(9)] = x4f;
    e[CINDEX(10)] = x5f;
}

///////////////////////////////////////////////////////////////////////////
///////////////////  Define the phase linkages function ///////////////////
///////////////////////////////////////////////////////////////////////////

void linkages(adouble* linkages, adouble* xad, Workspace* workspace)
{
    // No linkages as this is a single phase problem
}

////////////////////////////////////////////////////////////////////////////
///////////////////  Define the main routine ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(void)
{

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare key structures ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    Alg algorithm;
    Sol solution;
    Prob problem;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem name  ////////////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.name = "Time-Optimal Trajectory Generation Problem";

    problem.outfilename = "traj.txt";

////////////////////////////////////////////////////////////////////////////
////////////  Define problem level constants & do level 1 setup ////////////
////////////////////////////////////////////////////////////////////////////

    problem.nphases = 1;
    problem.nlinkages = 0;

    psopt_level1_setup(problem);

/////////////////////////////////////////////////////////////////////////////
/////////   Define phase related information & do level 2 setup /////////////
/////////////////////////////////////////////////////////////////////////////

    problem.phases(1).nstates = 5;
    problem.phases(1).ncontrols = 2;
    problem.phases(1).nevents = 10;
    problem.phases(1).npath = 2;
    problem.phases(1).nodes = "[100]";

    psopt_level2_setup(problem, algorithm);

////////////////////////////////////////////////////////////////////////////
///////////////////  Declare DMatrix objects to store results //////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix x, u, t;
    DMatrix lambda, H;

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter problem bounds information //////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.phases(1).bounds.lower.states(1) = -100.0;
    problem.phases(1).bounds.lower.states(2) = -100.0;
    problem.phases(1).bounds.lower.states(3) = -1 * M_PI;
    problem.phases(1).bounds.lower.states(4) = -10 * M_PI;
    problem.phases(1).bounds.lower.states(5) = -10 * M_PI;

    problem.phases(1).bounds.upper.states(1) = 5.0;
    problem.phases(1).bounds.upper.states(2) = 5.0;
    problem.phases(1).bounds.upper.states(3) = 1 * M_PI;
    problem.phases(1).bounds.upper.states(4) = 10 * M_PI;       // 5 rps
    problem.phases(1).bounds.upper.states(5) = 10 * M_PI;

    problem.phases(1).bounds.lower.controls(1) = -0.5;
    problem.phases(1).bounds.lower.controls(2) = -0.5;

    problem.phases(1).bounds.upper.controls(1) = 0.5;
    problem.phases(1).bounds.upper.controls(2) = 0.5;

    problem.phases(1).bounds.lower.events(1) = 0.0;
    problem.phases(1).bounds.lower.events(2) = 0.0;
    problem.phases(1).bounds.lower.events(3) = 0.0;
    problem.phases(1).bounds.lower.events(4) = 0.0;
    problem.phases(1).bounds.lower.events(5) = 0.0;
    problem.phases(1).bounds.lower.events(6) = 2.0;     // final x
    problem.phases(1).bounds.lower.events(7) = 2.0;     // final y
    problem.phases(1).bounds.lower.events(8) = M_PI;     // final theta
    problem.phases(1).bounds.lower.events(9) = 0.0;     // final d/dt phi_1
    problem.phases(1).bounds.lower.events(10) = 0.0;    // final d/dt phi_2

    problem.phases(1).bounds.upper.events = problem.phases(1).bounds.lower.events;

    problem.phases(1).bounds.lower.path(1) = -M_PI / 2;
    problem.phases(1).bounds.upper.path(1) = M_PI / 2;
    problem.phases(1).bounds.lower.path(2) = 0.5;
    problem.phases(1).bounds.upper.path(2) = 100;

    problem.phases(1).bounds.lower.StartTime = 0.0;
    problem.phases(1).bounds.upper.StartTime = 0.0;

    problem.phases(1).bounds.lower.EndTime = 0.0;
    problem.phases(1).bounds.upper.EndTime = 100.0;

////////////////////////////////////////////////////////////////////////////
///////////////////  Register problem functions  ///////////////////////////
////////////////////////////////////////////////////////////////////////////

    problem.integrand_cost = &integrand_cost;
    problem.endpoint_cost = &endpoint_cost;
    problem.dae = &dae;
    problem.events = &events;
    problem.linkages = &linkages;

////////////////////////////////////////////////////////////////////////////
///////////////////  Define & register initial guess ///////////////////////
////////////////////////////////////////////////////////////////////////////

    int nnodes = (int) problem.phases(1).nodes(1);

    DMatrix x0(5, nnodes);

    x0(1, colon()) = linspace(0.0, 10.0, nnodes);
    x0(2, colon()) = linspace(0.0, 0.0, nnodes);
    x0(3, colon()) = linspace(0.0, 0, nnodes);
    x0(4, colon()) = linspace(0.0, 0.0, nnodes);
    x0(5, colon()) = linspace(0.0, 0.0, nnodes);
    //x0(6, colon()) = linspace(0.0, 0.0, nnodes);

    DMatrix u0(2, nnodes);
    u0(1, colon()) = linspace(1.0, -1.0, nnodes);
    u0(2, colon()) = linspace(1.0, -1.0, nnodes);

    problem.phases(1).guess.controls = u0;
    problem.phases(1).guess.states = x0;
    problem.phases(1).guess.time = linspace(0.0, 100.0, nnodes);

////////////////////////////////////////////////////////////////////////////
///////////////////  Enter algorithm options  //////////////////////////////
////////////////////////////////////////////////////////////////////////////

    algorithm.nlp_method = "IPOPT";
    algorithm.scaling = "automatic";
    algorithm.derivatives = "automatic";
    algorithm.nlp_iter_max = 5000;
    algorithm.nlp_tolerance = 1.e-6;

////////////////////////////////////////////////////////////////////////////
///////////////////  Now call PSOPT to solve the problem   /////////////////
////////////////////////////////////////////////////////////////////////////

#if 0
    psopt(solution, problem, algorithm);

    if (solution.error_flag) exit(0);

////////////////////////////////////////////////////////////////////////////
///////////  Extract relevant variables from solution structure   //////////
////////////////////////////////////////////////////////////////////////////

    x = solution.get_states_in_phase(1);
    u = solution.get_controls_in_phase(1);
    t = solution.get_time_in_phase(1);

////////////////////////////////////////////////////////////////////////////
///////////  Save solution data to files if desired ////////////////////////
////////////////////////////////////////////////////////////////////////////

    x.Save("x.dat");
    u.Save("u.dat");
    t.Save("t.dat");

////////////////////////////////////////////////////////////////////////////
///////////  Plot some results if desired (requires gnuplot) ///////////////
////////////////////////////////////////////////////////////////////////////

    DMatrix pos_x = x(1, colon());
    DMatrix pos_y = x(2, colon());

    DMatrix alpha = colon(0.0, pi/20, 2*pi);

    DMatrix xObs1 = sqrt(0.5)*cos(alpha) + 1;
    DMatrix yObs1 = sqrt(0.5)*sin(alpha) + 1;

    plot(pos_x, pos_y, xObs1, yObs1, problem.name + ": x-y trajectory", "x", "y", "pos");

    plot(t, x, problem.name + ": states", "time (s)", "states", "x1 x2 x3 x4 x5 x6");

    plot(t, u, problem.name + ": control", "time (s)", "control", "u1 u2");

    plot(pos_x, pos_y, xObs1, yObs1, problem.name + ": x-y trajectory", "x", "y", "pos obs1", "pdf", "lts_xy.pdf");

    plot(t, x, problem.name + ": states", "time (s)", "states", "x1 x2 x3 x4 x5 x6", "pdf", "lts_states.pdf");

    plot(t, u, problem.name + ": control", "time (s)", "control", "u1 u2", "pdf", "lts_control.pdf");

#endif
    ////////////////////////////////////////////////////////////////////////////////
      // GLFW の初期化 (GLFW)
      if (glfwInit() == GL_FALSE){
        // 初期化に失敗した処理
        std::cerr << "Can't initialize GLFW" << std::endl;
        return 1;
      }


      ////////////////////////////////////////////////////////////////////////////////
      // 終了時の処理登録 (GLFW)
      atexit(glfwTerminate);


      ////////////////////////////////////////////////////////////////////////////////
      // ウィンドウを作成 (GLFW)
      GLFWwindow * const window(glfwCreateWindow(/* int           width   = */ 640,
                                                 /* int           height  = */ 480,
                                                 /* const char  * title   = */ "Hello!",
                                                 /* GLFWmonitor * monitor = */ NULL,
                                                 /* GLFWwindow  * share   = */ NULL));
      if (window == NULL){
        // ウィンドウ作成に失敗した処理
        std::cerr << "Can't create GLFW window." << std::endl;
        return 1;
      }

      // 作成したウィンドウを処理対象とする (GLFW)
      glfwMakeContextCurrent(/* GLFWwindow *  window = */ window);

      // 背景色 (OpenGL)
      glClearColor(/* GLfloat red   = */ 0.2f,
                   /* GLfloat green = */ 0.2f,
                   /* GLfloat blue  = */ 0.2f,
                   /* GLfloat alpha = */ 0.0f);


      ////////////////////////////////////////////////////////////////////////////////
      // ループ処理

      // ウィンドウが開いている間繰り返す
      while (glfwWindowShouldClose(/* GLFWwindow * window = */ window) == GL_FALSE){
        // ウィンドウを消去 (GLFW)
        glClear(/* GLbitfield mask = */ GL_COLOR_BUFFER_BIT);

        // 描画処理

        // カラーバッファ入れ替え <= ダブルバッファリング (GLFW)
        glfwSwapBuffers(window);

        // イベント待ち (GLFW)
        glfwWaitEvents();
      }
}

////////////////////////////////////////////////////////////////////////////
///////////////////////      END OF FILE     ///////////////////////////////
////////////////////////////////////////////////////////////////////////////

