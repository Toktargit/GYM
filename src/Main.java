public class Main {
    final static double ratio_T = 1;//разница температур для каллорического идеального воздуха
    final static double R = 287;
    final static double gamma = 1.4;
    final static double aFlow = 340.28;//speed of sound m/s
    final static double Cv = R / (gamma - 1);
    final static double Cp = gamma * Cv;
    final static double CourantNumber = 0.8; //0.5-0.8 btwn
    final static double LHORI = 10e-5;//plate length m
    final static double PFlow = 101325;//pressure pasca
    final static double TFlow = 288.16;//value of temperature see level
    final static double RoFlow = PFlow/(R*TFlow);
    final static double mah = 4;//mah number
    final static double miu0 = 1.7894 * 10e-5;
    final static double T0 = 288.16;
    final static double UFlow = mah * aFlow;

    final static double Re = RoFlow*UFlow*aFlow*LHORI / dynvis();
    final static double sigma = (5 * LHORI) / Math.sqrt(Re);
    final static double LVERT = 5*sigma;
    final static int Imax = 71;
    final static int Jmax = 71;
    static double[][] Miuflow = new double[Imax][Jmax];
    static double[][] k = new double[Imax][Jmax];
    final static double Pr = 0.71;//Prndtle number
    final static double dx = LHORI / (Imax - 1);
    final static double dy = LVERT / (Jmax - 1);
    static boolean converged = false;


    static double[][] tauxx = new double[Imax][Jmax];
    static double[][] tauxy = new double[Imax][Jmax];
    static double[][] tauyy = new double[Imax][Jmax];
    static double[][] qx = new double[Imax][Jmax];
    static double[][] qy = new double[Imax][Jmax];

    public static double max(double[][] mat){
        double max = 0.0;
        for (double[] doubles : mat) {
            max = doubles[0];
            for (int j = 1; j < doubles.length; j++) {
                if (doubles[j] > max) {
                    max = doubles[j];
                }
            }
        }
        return max;
    }

    static double[][] dynvis(double[][] T){
        for (int i = 0; i < Imax; i++) {
            for (int j = 0; j < Jmax; j++) {
                Miuflow[i][j] = Main.miu0 * (Math.pow(T[i][j], 3.0 / 2.0) / Math.pow(Main.T0, 3.0 / 2.0)) * ((Main.T0 + 110) / (T[i][j] + 110));
            }
        }
        return Miuflow;
    }
    static double dynvis(){
        return Main.miu0 * (Math.pow(Main.TFlow, 3.0 / 2.0) / Math.pow(Main.T0, 3.0 / 2.0)) * ((Main.T0 + 110) / (Main.TFlow + 110));
    }
    static double[][] thermc(double[][] miu){
        for (int i = 0; i < Imax; i++) {
            for (int j = 0; j < Jmax; j++) {
                k[i][j] = miu[i][j] * Main.Cp / Main.Pr;
            }
        }
        return k;
    }


    public static double[][] qy(double[][] T, double[][] k, double dy, String call_case){
        double[][] dT_dy = new double[Jmax][Imax];

        if ((call_case.equals("Predict_F"))){
            for (int i = 0; i < Imax-1; i++) {
                for (int j = 1; j < Jmax-1; j++) {
                    dT_dy[j][i] = (T[j][i] - T[j-1][i])/dy;
                }
            }
            for (int k1 = 0; k1 < Jmax-1; k1++) {
                dT_dy[0][k1] = (T[1][k1] - T[0][k1])/dy; // Forward at i = 1
            }
        }

        else if(call_case.equals("Correct_F")) {
            for (int i = 0; i < Imax-1; i++) {
                for (int j = 0; j < Jmax-2; j++) {
                    dT_dy[j][i] = (T[j + 1][i] - T[j][i])/dy;
                }
            }
            for (int k1 = 0; k1 < Jmax; k1++) {
                dT_dy[Jmax-1][k1] = (T[Jmax-1][k1] - T[Jmax - 2][k1]) / dy;
            }
        }

        else {
            throw new RuntimeException("Undefined call case.");
        }

        for (int i = 0; i < Imax-1; i++) {
            for (int j = 0; j < Jmax-1; j++) {
                qy[i][j] = - k[i][j]*dT_dy[i][j];
            }
        }
        return qy;
    }
    public static double[][] q_x (double[][] T, double[][] k, double dx, String call_case){
        double[][] dT_dx = new double[Jmax][Imax];

        if ((call_case.equals("Predict_E"))){
            for (int i = 1; i < Imax-1; i++) {
                for (int j = 0; j < Jmax-1; j++) {
                    dT_dx[j][i] = (T[j][i] - T[j][i - 1])/dx;
                }
            }
            for (int k1 = 0; k1 < Jmax-1; k1++) {
                dT_dx[k1][0] = (T[k1][1] - T[k1][0])/dx; // Forward at i = 1
            }
        }

        else if(call_case.equals("Correct_E")) {
            for (int i = 0; i < Imax-2; i++) {
                for (int j = 0; j < Jmax-1; j++) {
                    dT_dx[j][i] = (T[j][i + 1] - T[j][i])/dx;
                }
            }
            for (int k1 = 0; k1 < Jmax; k1++) {
                dT_dx[k1][Imax - 1] = (T[k1][Imax - 1] - T[k1][Imax - 2]) / dx;
            }
        }
        else {
            throw new RuntimeException("Undefined call case.");
        }

        for (int i = 0; i < Imax-1; i++) {
            for (int j = 0; j < Jmax-1; j++) {
                qx[i][j] = - k[i][j]*dT_dx[i][j];
            }
        }

        return qx;
    }


    static double[][] tauxx(double[][] u,double[][] v,double[][] lambda, double[][] miu, double dx,double dy, String call_case){

        double[][] du_dx = new double[Jmax][Imax];
        double[][] dv_dy = new double[Jmax][Imax];

        if ((call_case.equals("Predict_E"))){
            for (int i = 1; i < Imax-1; i++) {
                for (int j = 0; j < Jmax-1; j++) {
                    du_dx[j][i] = (u[j][i] - u[j][i-1])/dx;
                }
            }
            for (int k = 0; k < Jmax-1; k++) {
                du_dx[k][0] = (u[k][1] - u[k][0])/dx; // Forward at i = 1
            }
        }
        else if(call_case.equals("Correct_E"))
        {
            for (int i = 1; i < Imax-1; i++) {
                for (int j = 0; j < Jmax-1; j++) {
                    du_dx[j][i] = (u[j][i] - u[j][i-1])/dx;
                }
            }
            for (int k = 0; k < Jmax; k++) {
                du_dx[k][Imax-1] = u[k][Imax-1] /dx;
            }
        }
        else
        {
            throw new RuntimeException("Undefined call case.");
        }
        for (int i = 0; i < Imax-1; i++) {
            for (int j = 1; j < Jmax-2; j++) {
                dv_dy[j][i] = (v[j+1][i] - v[j-1][i]) / (2*dy);
            }
        }
        for (int k = 0; k < Imax; k++) {
            dv_dy[0][k] = (v[1][k] - v[0][k])/dy;
            dv_dy[Jmax-1][k] = (v[Jmax-1][k] - v[Jmax-2][k])/dy;
        }
        for (int i = 0; i < Imax; i++) {
            for (int j = 0; j < Jmax; j++) {
                tauxx[j][i] = lambda[j][i] * (du_dx[j][i] + dv_dy[j][i]) + 2 * miu[j][i] * du_dx[j][i];
            }
        }
        return tauxx;
    }
    static double[][] tauxy(double[][] u,double[][] v, double[][] miu, double dx,double dy, String call_case) {

        double[][] du_dy = new double[Jmax][Imax];
        double[][] dv_dx = new double[Jmax][Imax];
        if ((call_case.equals("Predict_E")) || (call_case.equals("Correct_E"))) {
            for (int i = 0; i < Imax - 1; i++) {
                for (int j = 1; j < Jmax - 2; j++) {
                    du_dy[j][i] = (u[j + 1][i] - u[j - 1][i]) / (2 * dy);
                }
            }
            for (int k = 0; k < Imax-1; k++) {
                du_dy[0][k] = (u[1][k] - u[0][k]) / dy;
                du_dy[Jmax-1][k] = (u[Jmax-1][k] - u[Jmax - 2][k]) / dy;
            }
            if ((call_case.equals("Predict_E"))) {
                for (int i = 1; i < Imax - 1; i++) {
                    for (int j = 0; j < Jmax - 1; j++) {
                        dv_dx[j][i] = (v[j][i] - v[j][i - 1]) / dx;
                    }
                }
                for (int k = 0; k < Jmax - 1; k++) {
                    dv_dx[k][0] = (v[k][1] - v[k][0]) / dx; // Forward at i = 1
                }
            } else {
                for (int i = 0; i < Imax - 2; i++) {
                    for (int j = 0; j < Jmax - 1; j++) {
                        dv_dx[j][i] = (v[j][i + 1] - v[j][i]) / dx;
                    }
                }
                for (int i = 0; i < Jmax - 1; i++) {
                    dv_dx[i][Imax - 1] = (v[i][Imax - 1] - v[i][Imax - 2]) / dx;
                }
            }
        }
        else if ((call_case.equals("Predict_F")) || (call_case.equals("Corrector_F"))) {
            for (int i = 1; i < Imax - 2; i++) {
                for (int j = 0; j < Jmax - 1; j++) {
                    dv_dx[j][i] = (v[j][i + 1] - v[j][i - 1]) / (2 * dx);
                }
            }
            for (int i = 0; i < Jmax; i++) {
                dv_dx[i][0] = (v[i][1] - v[i][0]) / dx;
                dv_dx[i][Imax - 1] = (v[i][Imax - 1] - v[i][Imax - 2]) / dx;
            }
            if ((call_case.equals("Predict_F"))) {
                for (int i = 0; i < Imax - 1; i++) {
                    for (int j = 1; j < Jmax - 1; j++) {
                        du_dy[j][i] = (u[j][i] - u[j - 1][i]) / dy;
                    }
                }
                for (int i = 0; i < Imax - 1; i++) {
                    du_dy[0][i] = (u[1][i] - u[0][i]) / dy;
                }
            } else {
                for (int i = 0; i < Imax - 1; i++) {
                    for (int j = 0; j < Jmax - 2; j++) {
                        du_dy[j][i] = (u[j + 1][i] - u[j][i]) / dy;
                    }
                }
                for (int i = 0; i < Jmax - 1; i++) {
                    du_dy[Jmax-1][i] = (u[Jmax-1][i] - u[Jmax - 2][i]) / dy;
                }
            }
        }
        else {
            throw new RuntimeException("Undefined call case.");
        }
        for (int i = 0; i < Imax-1; i++) {
            for (int j = 0; j < Jmax-1; j++) {
                tauxy[i][j] = miu[i][j] * (du_dy[i][j] + dv_dx[i][j]);
            }
        }

        return tauxy;
    }
    static double[][] tauyy(double[][] u,double[][] v,double[][] lambda, double[][] miu, double dx,double dy, String call_case){

        double[][] du_dx = new double[Jmax][Imax];
        double[][] dv_dy = new double[Jmax][Imax];

        if ((call_case.equals("Predict_F"))){
            for (int i = 0; i < Imax-1; i++) {
                for (int j = 1; j < Jmax-1; j++) {
                    dv_dy[j][i] = (v[j][i] - v[j-1][i])/dy;
                }
            }
            for (int k = 0; k < Jmax-1; k++) {
                dv_dy[0][k] = (v[1][k] - v[0][k])/dy; // Forward at i = 1
            }
        }
        else if(call_case.equals("Correct_F")) {
            for (int i = 0; i < Imax-1; i++) {
                for (int j = 0; j < Jmax-2; j++) {
                    dv_dy[j][i] = (v[j+1][i] - v[j][i])/dy;
                }
            }
            for (int k = 0; k < Imax-1; k++) {
                dv_dy[Jmax - 1][k] = (v[Jmax - 1][k] - v[Jmax - 2][k])/dy;
            }
        }
        else {
            throw new RuntimeException("Undefined call case.");
        }
        for (int i = 1; i < Imax-2; i++) {
            for (int j = 0; j < Jmax-1; j++) {
                du_dx[j][i] = (u[j][i + 1] - u[j][i - 1]) / (2*dx);
            }
        }
        for (int k = 0; k < Jmax; k++) {
            du_dx[k][0] = (u[k][1] - u[k][0])/dx;
            du_dx[k][Imax - 1] = (u[k][Imax - 1] - u[k][Imax - 2])/dx;
        }
        for (int i = 0; i < Imax; i++) {
            for (int j = 0; j < Jmax; j++) {
                tauyy[j][i] = lambda[j][i] * (du_dx[j][i] + dv_dy[j][i]) + 2 * miu[j][i] * dv_dy[j][i];
            }
        }

        return tauyy;
    }


    public static void intermediateE (double[][] Ro, double[][] u,double [][]  p,double[][] v,double[][] T,double[][] miu,double[][] lambda,double [][] k,double Cv,double dx,double dy, String call_case){
        tauxy = tauxy(u, v, miu, dx, dy, call_case);
        tauxx = tauxx(u, v, lambda, miu, dx, dy, call_case);
        for (int i = 0; i < Imax-1; i++) {
            for (int j = 0; j < Jmax-1; j++) {
                E1[i][j] = Ro[i][j] * u[i][j];
                //tauxy
                E2[i][j] = Math.pow(Ro[i][j] * u[i][j],2) + p[i][j] - tauxx[i][j];
                //tauxx
                E3[i][j] = Ro[i][j] * u[i][j] * v[i][j] - tauxy[i][j];
                qx = q_x(T, k, dx, call_case);
                E5[i][j] = (Ro[i][j] * (Cv*T[i][j] + (Math.pow(u[i][j],2) + Math.pow(v[i][j],2))/2) + p[i][j])* u[i][j] - u[i][j]*tauxx[i][j] - v[i][j]*tauxy[i][j] + qx[i][j];
            }
        }
    }
    public static void intermediateF (double[][] Ro, double[][] u,double [][]  p,double[][] v,double[][] T,double[][] miu,double[][] lambda,double [][] k,double Cv,double dx,double dy, String call_case){
        tauxy = tauxy(u, v, miu, dx, dy, call_case);
        tauyy = tauyy(u, v, lambda, miu, dx, dy, call_case);
        for (int i = 0; i < Imax - 1; i++) {
            for (int j = 0; j < Jmax - 1; j++) {
                F1[i][j] = Ro[i][j] * v[i][j];
                //tauxy = tauxy(u, v, miu, dx, dy, "call_case");
                F2[i][j] = Ro[i][j] * u[i][j] * v[i][j] - tauxy[i][j];
               // tauyy = tauyy(u, v, lambda, miu, dx, dy, "call_case");
                F3[i][j] = Math.pow(Ro[i][j] * v[i][j],2) +p[i][j] - tauyy[i][j];
                qy = qy(T, k, dy, call_case);
                F5[i][j] = (Ro[i][j] * (Cv*T[i][j] + (Math.pow(u[i][j],2) + Math.pow(v[i][j],2))/2) + p[i][j]) * v[i][j] - u[i][j]*tauxy[i][j] - v[i][j]*tauyy[i][j] + qy[i][j];
            }
        }
    }


    static double[][] E1 = new double[Imax][Jmax];
    static double[][] E2 = new double[Imax][Jmax];
    static double[][] E3 = new double[Imax][Jmax];
    static double[][] E5 = new double[Imax][Jmax];

    static double[][] F1 = new double[Imax][Jmax];
    static double[][] F2 = new double[Imax][Jmax];
    static double[][] F3 = new double[Imax][Jmax];
    static double[][] F5 = new double[Imax][Jmax];

    static double[][] Ro = new double[Imax][Jmax];
    static double[][] u = new double[Imax][Jmax];
    static double[][] v = new double[Imax][Jmax];
    static double[][] p = new double[Imax][Jmax];
    static double[][] T = new double[Imax][Jmax];

    public static void BC(double[][] Ro_p, double[][] u_p, double[][] v_p, double[][] p_p, double[][] T_p, double RoFlow, double UFlow, double PFlow, double TFlow) {
        //leading edge
        T[0][0] = TFlow;
        p[0][0] = PFlow;
        Ro[0][0] = RoFlow;
        //Inflow / Upper boundary (not leading edge) --> y component of velocity (v) is assumed equal to zero

//% v(2:JMAX,1) = 0 is implicit in the definition of v
        for (int j = 1; j < Jmax-1; j++) {
                u[j][0] = UFlow;
                p[j][0] = PFlow;
                T[j][0] = TFlow;
                Ro[j][0] = RoFlow;
            }
       // For the upper boundary:
        for (int i = 0; i < Imax-1; i++) {
            u[Jmax-1][i] = UFlow;
            p[Jmax-1][i] = PFlow;
            T[Jmax-1][i] = TFlow;
            Ro[Jmax-1][i] = RoFlow;
        }
//% Case 4: Outflow (not surface or JMAX) --> Velocity, pressure and
           //     % temperature are calculated based on an extrapolation from the two adjacent interior points.
// Density is computed from the equation of state:
        for (int i = 1; i < Jmax-2; i++) {
            u[i][Imax-1] = 2*u[i][Imax-2] - u[i][Imax-3];
            v[i][Imax-1] = 2*v[i][Imax-2] - v[i][Imax-3];
            p[i][Imax-1] = 2*p[i][Imax-2] - p[i][Imax-3];
            T[i][Imax-1] = 2*T[i][Imax-2] - T[i][Imax-3];
            Ro[i][Imax-1] = p[i][Imax-1] / (R*T[i][Imax-1]);
        }
        // Surface (not leading edge) --> No-slip condition is specified on velocity. u[0][i]=0 v[0][i] = 0
        for (int i = 1; i < Imax-1; i++) {
            T[0][i] = ratio_T * TFlow;
            p[0][i] = 2*p[1][i] - p[2][i];
            Ro[0][i] = p[0][i] / (R*T[0][i]);
        }


    }

    public static boolean conv(double[][] Ro_old, double[][] Ro){
        double[][] result = new double[Imax][Jmax];
        for (int i = 0; i < Imax; i++) {
            for (int j = 0; j < Jmax; j++) {
                result[i][j] = Math.abs(Ro_old[i][j] - Ro[i][j]);
                if(max(result) < 1e-8){
                    converged = true;
                }
                else{
                    converged = false;
                }
                System.out.print(Ro_old[i][j] + "   "  );
            }
            System.out.println();
        }

        return converged;
    }

    public static void main(String[] args) {

        int Maxiter = 4000, t = 1; //time step counter
        double time = 0;//physical time

        //double EFlow  = Cv * TFlow; специфичная энергия, ненужен
        // kFlow = thermc(Pr, cp, miuFlow);специфичная теплопроводимость, ненужен

        double[] x = new double[Imax];
        for (int i = 0; i < Imax; i++) {
            if (x[i] == LHORI) break;
            x[i] += dx;
        }

        double[] y = new double[Jmax];
        for (int i = 0; i < Jmax; i++) {
            if (y[i] == LVERT) break;
            y[i] += dy;
        }
        //переменные для нахождение ответов уровнение U E F вот этих скорее всего надо добавить в массив W

        double[][] miu = dynvis(T);
        double[][] lambda = new double[Imax][Jmax];
        double[][] k = thermc(miu);

        for (int i = 0; i < Imax; i++) {
            for (int j = 0; j < Jmax; j++) {
                p[i][j] = 1*PFlow;
                Ro[i][j] = 1*RoFlow;
                T[i][j] = 1*TFlow;
                u[i][j] = 1*UFlow;

                lambda[i][j] = -2.0/3.0 * miu[i][j];

            }
        }
        for (int i = 0; i < Jmax; i++) {
            T[0][i] = ratio_T*TFlow;
            u[0][i] = 0;
        }

        double[][] U1_p = new double[Imax][Jmax];
        double[][] U2_p = new double[Imax][Jmax];
        double[][] U3_p = new double[Imax][Jmax];
        double[][] U5_p = new double[Imax][Jmax];
        double[][] Ro_p = new double[Imax][Jmax];
        double[][] u_p = new double[Imax][Jmax];
        double[][] v_p = new double[Imax][Jmax];
        double[][] T_p = new double[Imax][Jmax];
        double[][] p_p = new double[Imax][Jmax];

        double vSPEED;

        double[][][] U = new double[4][Imax][Jmax];
        double[][][] E = new double[4][Imax][Jmax];
        double[][][] F = new double[4][Imax][Jmax];

        //start

        double[][] Ro_old =new double[Imax][Jmax];//= Ro;
        for (int i = 0; i < Imax; i++) {
            for (int j = 0; j < Jmax; j++) {
                Ro_old[i][j] = Ro[i][j];
            }
        }
        while(!converged && t <= Maxiter) {
            double[][] temp = new double[Imax][Jmax];
            double max = temp[0][0];
            for (int i = 1; i < Imax - 2; i++) {
                for (int j = 1; j < Jmax - 2; j++) {
                    temp[i][j] = ((4.0 / 3.0 * Miuflow[i][j] * (gamma * Miuflow[i][j] / Pr)) / Ro[i][j]);
                    if (temp[i][j] > max) {
                        max = temp[i][j];
                    }
                }
            }
            vSPEED = max;
//нахождение в с черточкой

            double[][] DTcfl = new double[Imax][Jmax];

            double min = 1e-6;
            for (int i = 1; i < Imax-2; i++) {// со второго до предпоследнего
                for (int j = 1; j < Jmax-2; j++) {
                    DTcfl[i][j] = 1 / Math.abs(u[i][j])/ dx +  Math.abs(v[i][j]) / dy
                            + Math.sqrt( gamma * R * T[i][j])
                            * Math.sqrt((1 / (dx * dx)) + (1 / (dy * dy)))
                            + 2 * vSPEED * Math.sqrt((1 / (dx * dx)) + (1 / (dy * dy)));
                    if (DTcfl[i][j] < min) {
                        min = DTcfl[i][j];
                    }
                }
            }
            double[] dt = new double[Maxiter];
            dt[t] = CourantNumber * min;

            //maccormac
            //решение юшки по времени
            double[][] U1 = new double[Imax][Jmax];
            double[][] U2 = new double[Imax][Jmax];
            double[][] U3 = new double[Imax][Jmax];
            double[][] U5 = new double[Imax][Jmax];

            for (int i = 0; i < Imax; i++) {
                for (int j = 0; j < Jmax; j++) {
                    U1[i][j] = Ro[i][j]; // From continuity
                    U2[i][j] = Ro[i][j] * u[i][j]; // From momentum-x умножаем поэлементно
                    U3[i][j] = Ro[i][j] * v[i][j]; // From momentum-y
                    U5[i][j] = Ro[i][j] * (Cv * T[i][j] + (Math.pow(u[i][j], 2) + Math.pow(v[i][j], 2))); // From energy
                }
            }

            //Predictor
                intermediateE(Ro, u, p, v, T, miu, lambda, k, Cv, dx, dy, "Predict_E");
                intermediateF(Ro, u, p, v, T, miu, lambda, k, Cv, dx, dy,"Predict_F");

            for (int i = 2; i < Imax-1; i++) {
                for (int j = 2; j < Jmax-1; j++) {
                    U1_p[j][i] =  U1_p[j][i] - dt[t] * ((E1[j][i+1] - E1[j][i])/dx + (F1[j+1][i] - F1[j][i])/dy);
                    U2_p[j][i] =  U2_p[j][i] - dt[t] * ((E2[j][i+1] - E2[j][i])/dx + (F2[j+1][i] - F2[j][i])/dy);
                    U3_p[j][i] =  U3_p[j][i] - dt[t] * ((E3[j][i+1] - E3[j][i])/dx + (F3[j+1][i] - F3[j][i])/dy);
                    U5_p[j][i] =  U5_p[j][i] - dt[t] * ((E5[j][i+1] - E5[j][i])/dx + (F5[j+1][i] - F5[j][i])/dy);
                }
            }
            //[1 до JMAX-1][IMAX-1]
            //u_p[1 до JMAX-1][IMAX-1]
            //V_P[1 до JMAX-1][IMAX-1]
            //T_p[1 до JMAX-1][IMAX-1]

//            double[][][] W = new double[4][Imax][Jmax];
//            for (int i = 1; i < Imax-2; i++) {
//                for (int j = 1; j < Jmax-2; j++) {
//                    W[0][i][j] = U1[i][j];//ro
//                    W[1][i][j] = U2[i][j]/U1[i][j];//u
//                    W[2][i][j] = U3[i][j]/U1[i][j];//v
//                    W[3][i][j] = (U5[i][j]/U1[i][j] - (Math.pow((U2[i][j]/U1[i][j]),2) + Math.pow((U3[i][j]/U1[i][j]), 2))/2)/Cv;//T
//                }
//            }
            double[][][] W = new double[4][Imax][Jmax];
            for (int i = 1; i < Imax-2; i++) {
                for (int j = 1; j < Jmax - 2; j++) {
                    Ro_p[i][j] = U1[i][j];                //ro
                    u_p[i][j] = U2[i][j] / U1[i][j];      //u
                    v_p[i][j] = U3[i][j] / U1[i][j];      //v
                    T_p[i][j] = (U5[i][j] / U1[i][j] - (Math.pow((U2[i][j] / U1[i][j]), 2) + Math.pow((U3[i][j] / U1[i][j]), 2)) / 2) / Cv;//T
                }
            }

            for (int i = 1; i < Imax-2; i++) {
                for (int j = 1; j < Jmax-2; j++) {
                    p_p[j][i] = Ro_p[j][i] * R * T_p[j][i];
                }
            }

            //static [Ro_p, u_p, v_p, p_p, T_p]
             BC(Ro_p, u_p, v_p, p_p, T_p, RoFlow , UFlow, PFlow, TFlow);

            converged = conv(Ro_old, Ro);
            Ro = Ro_old;
            time = time + dt[t];
            t++;
            System.out.println(time);
        }
    }
}
