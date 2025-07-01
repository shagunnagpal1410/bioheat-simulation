#include<iostream>
#include<vector>
#include<map>
#include<cmath>
#include<Eigen/Sparse>
#include<Eigen/Dense>
#include<Eigen/SparseLU>
#include<fstream>
#include<sstream>
#include<string>
using namespace std;
using namespace Eigen;
typedef Triplet<double> Tri;
vector<vector<double>> readCSV(const string& filename) {
    vector<vector<double>> data;
    ifstream file(filename);
    string line;
    getline(file,line);
    while (getline(file,line)) {
        stringstream ss(line);
        string val;
        vector<double> row;
        while (getline(ss,val,',')) {
            row.push_back(stod(val));
        }
        if (!row.empty()) {
            data.push_back(row);
        }
    }
    return data;
}
struct point{
    double x,y,z;
    int voxel; 
    point(double i, double j, double k) {
        x=i, y=j, z=k;
        voxel=0;
    }
    bool operator<(const point& other) const {
        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
        if (z != other.z) return z < other.z;
        return voxel < other.voxel;
    }
};
int provide_voxel(int x_number, int y_number, int z_number, int voxels_inrow, int voxels_incolumn) {
    return z_number*voxels_incolumn*voxels_inrow+y_number*voxels_inrow+x_number;
}
int provide_znumber(int voxel, int voxels_inrow, int voxels_incolumn) {
    return int(voxel/(voxels_inrow*voxels_incolumn));
}
int provide_ynumber(int voxel, int voxels_inrow, int voxels_incolumn) {
    return (voxel-provide_znumber(voxel,voxels_inrow, voxels_incolumn)*voxels_inrow*voxels_incolumn)/voxels_inrow;
}
int provide_xnumber(int voxel, int voxels_inrow, int voxels_incolumn) {
    return (voxel-provide_znumber(voxel,voxels_inrow, voxels_incolumn)*voxels_inrow*voxels_incolumn-provide_ynumber(voxel,voxels_inrow, voxels_incolumn)*voxels_inrow);
}
double find_Qm(point p1, double x1, double y1, double z1, double x2, double y2, double z2, double x3, double y3,double z3, double a1, double a2, double a3) {
    double numerator =0;
    double sigma=0.003;
    double denominator = 0;
    double x=p1.x, y=p1.y, z=p1.z;
    numerator+=a1*exp(-1*((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1))/(2*sigma*sigma));
    numerator+=a2*exp(-1*((x-x2)*(x-x2)+(y-y2)*(y-y2)+(z-z2)*(z-z2))/(2*sigma*sigma));
    numerator+=a3*exp(-1*((x-x3)*(x-x3)+(y-y3)*(y-y3)+(z-z3)*(z-z3))/(2*sigma*sigma));
    denominator+=(a1+a2+a3)*(pow(sigma*pow(2*M_PI,0.5),3));
    return ((numerator*16.6) / denominator);
}
double find_W(point Ni, point p) {
    double value=pow(Ni.x-p.x,2)+pow(Ni.y-p.y,2)+pow(Ni.z-p.z,2);
    double radius=10/1000.0;
    if (pow(value,0.5)/radius>1) {
        return 0.0;
    }
    else {
        return exp(-1*6.0*value/(radius*radius));
    }
}
int main() {
    double rhob=1060, Cb=4180, Ta=37.0;
    ofstream fout("6june.csv");
    fout<<"X"<<","<<"Y"<<","<<"Z"<<","<<"Temperature"<<"\n";
    double rhot=1030, Ct=3600, Kt=0.497;
    int T;
    double d;
    cout<<"Enter maximum time: ";
    cin>>T;
    cout<<"Enter time step: ";
    cin>>d;
    cout<<endl;
    int n=int(T/d);
    //asking for needle data
    double x1,y1,z1,x2,y2,z2,x3,y3,z3;
    cout<<"Enter x-coordinate of needle 1: ";
    cin>>x1;
    cout<<"Enter y-coordinate of needle 1: ";
    cin>>y1;
    cout<<"Enter z-coordinate of needle 1: ";
    cin>>z1;
    cout<<"Enter x-coordinate of needle 2: ";
    cin>>x2;
    cout<<"Enter y-coordinate of needle 2: ";
    cin>>y2;
    cout<<"Enter z-coordinate of needle 2: ";
    cin>>z2;
    cout<<"Enter x-coordinate of needle 3: ";
    cin>>x3;
    cout<<"Enter y-coordinate of needle 3: ";
    cin>>y3;
    cout<<"Enter z-coordinate of needle 3: ";
    cin>>z3;
    cout<<endl;
    double alpha1,alpha2, alpha3;
    cout<<"Enter alpha for needle 1: ";
    cin>>alpha1;
    cout<<"Enter alpha for needle 2: ";
    cin>>alpha2;
    cout<<"Enter alpha for needle 3: ";
    cin>>alpha3;
    cout<<endl;
    vector<vector<double>> volume_points=readCSV("volume2_points.csv");
    vector<vector<double>> tumor_points=readCSV("volume3_points.csv");
    for (int i=0; i<tumor_points.size(); i++) {
        volume_points.push_back({tumor_points[i][0], tumor_points[i][1], tumor_points[i][2]});
    }
    tumor_points.clear();
    vector<vector<double>> surface_points=readCSV("boundary_nodes.csv");
    //we will be doing origin shifting for more clarity
    for (int i=0; i<volume_points.size(); i++) {
        volume_points[i][0]+=104;
        volume_points[i][1]-=70;
        volume_points[i][2]+=118;
        volume_points[i][0]/=1000.0;
        volume_points[i][1]/=1000.0;
        volume_points[i][2]/=1000.0;
    }
    for (int i=0; i<surface_points.size(); i++) {
        surface_points[i][0]+=104;
        surface_points[i][1]-=70;
        surface_points[i][2]+=118;
        surface_points[i][0]/=1000.0;
        surface_points[i][1]/=1000.0;
        surface_points[i][2]/=1000.0;
    }
    //dimensions of cube 98 X 100 X 93
    double radius=10.0/1000.0;
    int voxels_inrow=int((98.0/1000.0)/radius)+1;
    int voxels_incolumn=int((100.0/1000.0)/radius)+1;
    int voxels_inz=int((93.0/1000.0)/radius)+1;
    vector<point> points;
    for (int i=0; i<volume_points.size(); i++) {
        point p1=point(volume_points[i][0],volume_points[i][1], volume_points[i][2]);
        p1.voxel=provide_voxel(int(volume_points[i][0]/radius), int(volume_points[i][1]/radius), int(volume_points[i][2]/radius), voxels_inrow, voxels_incolumn);
        points.push_back(p1);
    }
    int max_voxel=voxels_incolumn*voxels_inrow*voxels_inz;
    vector<vector<point>> points_insidevoxel(max_voxel);
    for (int i=0; i<points.size(); i++) {
        const point& p1=points[i];
        int voxel_num=p1.voxel;
        points_insidevoxel[voxel_num].push_back(p1);
    }
    vector<vector<int>> voxel_neighbours(max_voxel);
    for (int i=0; i<max_voxel; i++) {
        int x_number=provide_xnumber(i,voxels_inrow,voxels_incolumn);
        int y_number=provide_ynumber(i,voxels_inrow,voxels_incolumn);
        int z_number=provide_znumber(i,voxels_inrow,voxels_incolumn);
        for (int dx=-1; dx<=1; dx++) {
            for (int dy=-1; dy<=1; dy++) {
                for (int dz=-1; dz<=1; dz++) {
                    if (x_number+dx>=0 && x_number+dx<=voxels_inrow-1 && y_number+dy>=0 && y_number+dy<=voxels_incolumn-1 && z_number+dz>=0 && z_number+dz<=voxels_inz-1) {
                            voxel_neighbours[i].push_back(provide_voxel(x_number+dx, y_number+dy, z_number+dz, voxels_inrow, voxels_incolumn));
                        }
                }
            }
        }
    }
    map<point, double> checklist;
    for (int i=0; i<surface_points.size(); i++) {
        point p1={surface_points[i][0], surface_points[i][1], surface_points[i][2]};
        p1.voxel=provide_voxel(int(surface_points[i][0]/radius), int(surface_points[i][1]/radius), int(surface_points[i][2]/radius), voxels_inrow, voxels_incolumn);
        checklist[p1]=37.0;
    }
    int number=0;
    map<point, int> identity;
    for (int i=0; i<points.size(); i++) {
        point p1=points[i];
        if (checklist.find(p1)==checklist.end()) {
            identity[p1]=number;
            number+=1;
        }
    }
    int known_points=checklist.size();
    int unknown_points=identity.size();
    VectorXd u_prev=VectorXd :: Constant(unknown_points,37.0);
    vector<RowVectorXd> CT(points.size());
    cout<<points.size()<<endl;
    vector<vector<point>> unknown_neighbours_list(points.size());
    vector<vector<point>> known_neighbours_list(points.size());
    ofstream fout2("check_validity.csv");
    for (int p=0; p<points.size(); p++) {
        //iteration over points
        if(p%1000==0) cout<<"step- "<<p<<" done"<<endl;
            point p0=points[p];
            if (checklist.find(p0)!=checklist.end()) {
                fout2<<p<<","<<"known point"<<"\n";
                continue;
            }
            else {
                int total_neighbours=0;
                vector<point> known_neighbours;
                vector<point> unknown_neighbours;
                int voxel_number=p0.voxel;
                for (int i=0; i<voxel_neighbours[voxel_number].size(); i++) {
                    //neighbouring voxels starting
                    int neighbour_voxel=voxel_neighbours[voxel_number][i];
                    for (int j=0; j<points_insidevoxel[neighbour_voxel].size(); j++ ) {
                        //collecting neighbourhood points
                        point p2=points_insidevoxel[neighbour_voxel][j];
                        double distance=pow(p0.x-p2.x,2)+pow(p0.y-p2.y,2)+pow(p0.z-p2.z,2);
                        if (distance>0.0 && distance<pow(radius,2.0)) {
                            if(checklist.find(p2)!=checklist.end()) {
                                known_neighbours.push_back(p2);
                            }
                            else {
                                unknown_neighbours.push_back(p2);
                            }
                            total_neighbours+=1;
                        }
                        //neighbourhood points collected
                        if(total_neighbours>=200) break;
                    }
                    //neighbouring voxels ending
                    if(total_neighbours>=200) break;
                }
                total_neighbours=known_neighbours.size()+unknown_neighbours.size();
                MatrixXd M(total_neighbours,10);
                RowVectorXd L(10);
                MatrixXd W=MatrixXd :: Zero(total_neighbours, total_neighbours);
                number=0;
                for(int i=0; i<unknown_neighbours.size(); i++) {
                    point p3=unknown_neighbours[i];
                    double dx=p3.x-p0.x, dy=p3.y-p0.y, dz=p3.z-p0.z;
                    M(number,0)=1, M(number,1)=dx, M(number,2)=dy, M(number,3)=dz, M(number,4)=(dx*dx)/2, M(number,5)=(dy*dy)/2,
                    M(number,6)=(dz*dz)/2, M(number,7)=dx*dy, M(number,8)=dy*dz, M(number,9)=dz*dx;
                    W(number,number)=find_W(p3,p0);
                    number+=1;
                }
                for(int i=0; i<known_neighbours.size(); i++) {
                    point p3=known_neighbours[i];
                    double dx=p3.x-p0.x, dy=p3.y-p0.y, dz=p3.z-p0.z;
                    M(number,0)=1, M(number,1)=dx, M(number,2)=dy, M(number,3)=dz, M(number,4)=(dx*dx)/2, M(number,5)=(dy*dy)/2,
                    M(number,6)=(dz*dz)/2, M(number,7)=dx*dy, M(number,8)=dy*dz, M(number,9)=dz*dx;
                    W(number,number)=find_W(p3,p0);
                    number+=1;
                }
                L<<0,0,0,0,1,1,1,0,0,0;
                MatrixXd MTWM = M.transpose() * W * M;
                CT[p]=(L * MTWM.ldlt().solve(M.transpose() * W));
                known_neighbours_list[p]=known_neighbours;
                unknown_neighbours_list[p]=unknown_neighbours;
                fout2<<p<<","<<"unknown point"<<",";
                for (int i=0; i<CT[p].size(); i++) {
                    if (i!=CT[p].size()-1) {
                        fout2<<CT[p][i]<<",";
                    }
                    else {
                        fout2<<CT[p][i]<<"\n";
                    }
                }
                // Check size

            }
            //iterations over points ending
        }
        fout2.close();
    int store_time=INT_MAX;
    for (int t=1; t<=n; t++) {
    //starting of time loop
        cout<<"time step- "<<t<<" started!"<<endl;
        SparseMatrix<double> A(unknown_points,unknown_points);
        VectorXd rhs(unknown_points);
        vector<Tri> coefficients;
        int row_number;
        for (int p=0; p<points.size(); p++) {
        //iteration over points
            point p0=points[p];
            if (checklist.find(p0)!=checklist.end()) {
                continue;
            }
            else {
                //forming coefficient matrix
                row_number=identity[p0];
                double Qm=find_Qm(p0, x1,y1,z1,x2,y2,z2,x3,y3,z3,alpha1,alpha2,alpha3);
                if (u_prev(row_number)>=105) {
                    store_time=min(store_time,t);
                    Qm=0;
                }
                if (t>=store_time) {
                    Qm=0;
                }
                double wb=3.6*0.001;
                vector<point> known_neighbours=known_neighbours_list[p];
                vector<point> unknown_neighbours=unknown_neighbours_list[p];
                int total_neighbours=known_neighbours.size()+unknown_neighbours.size();
                double constant=(u_prev(identity[p0])*((-rhot*Ct/(Kt*d))))-(Qm/Kt)-(rhob*wb*Cb*Ta/Kt);
                for (int i=0; i<unknown_neighbours.size(); i++) {
                    point p3=unknown_neighbours[i];
                    coefficients.push_back(Tri(row_number, identity[p3], CT[p](i)));
                }
                for (int i=0; i<known_neighbours.size(); i++){
                    constant-=CT[p](unknown_neighbours.size()+i)*checklist[known_neighbours[i]];
                }
                coefficients.push_back(Tri(row_number, row_number, -(rhot*Ct/(Kt*d))-(rhob*wb*Cb/Kt)));
                rhs(row_number)=constant;
                //coeffecient matrix is done
            }
            //iterations over points ending
        }
            A.setFromTriplets(coefficients.begin(), coefficients.end());
            bool has_zero_row = false;
for (int i = 0; i < A.rows(); ++i) {
    bool nonzero_found = false;
    for (Eigen::SparseMatrix<double>::InnerIterator it(A, i); it; ++it) {
        if (std::abs(it.value()) > 1e-12) {
            nonzero_found = true;
            break;
        }
    }
    if (!nonzero_found) {
        cout << " Row " << i << " of matrix A is all zeros!" << endl;
        has_zero_row = true;
    }
}
if (has_zero_row) {
    cout << "Matrix A has one or more all-zero rows. SparseLU will fail." << endl;
}
            cout<<"Matrix Successfully formed!"<<endl;

            SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
            solver.compute(A);
            if (solver.info() !=Success) {
                cerr << "[ERROR] Decomposition failed!\n";
                return -1;
            }

            VectorXd u_next = solver.solve(rhs);
            if (solver.info() != Success) {
                cerr << "[ERROR] Solving failed!\n";
                return -1;
            }

            cout << "[SUCCESS] System solved! Max temp: " << u_next.maxCoeff() << "\n";
            u_prev = u_next;
    //time loop ending here
    }
    for (int i=0; i<points.size(); i++) {
        point p2=points[i];
        if (checklist.find(p2)!=checklist.end()) continue;
        int identity_number=identity[p2];
        checklist[p2]=u_prev(identity_number);
    }
    for (int i=0; i<points.size(); i++) {
        point p0=points[i];
        fout<<p0.x<<","<<p0.y<<","<<p0.z<<","<<checklist[p0]<<"\n";
    }
    fout.close();
    return 0;
}
