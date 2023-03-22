#include<iostream>
#include<algorithm>
#include<vector>
using namespace std;

typedef struct      //三维点的横纵列坐标
{                 
    double x;
    double y;
    double z;
}Point;

typedef struct      //三维点之间的边
{
    /* data */
    int left;
    int right;
    int count;      //用于控制边的计数，若为0，则删除此边
}Edge;

typedef struct{
    //圆心的坐标点
    double x_c;
    double y_c;
    double z_c;
    Point point[3];


}Triangle;




/*
三角剖分类
*/
class DelaunayTriangulation{
    public:
    vector<Point> point;                            //生成point数组
    vector<Triangle> triangle;                      //生成triangle数组
    vector<Edge> edge;
    DelaunayTriangulation(vector<Point> points);
    //构造初始大三角
    DelaunayTriangulation(){}                       //无参构造函数
    ~DelaunayTriangulation();                       //析构函数
    void Circumscribe(double &x_center,             //外接圆计算函数
                      double &y_center, 
                      double &r, 
                      Point v1, 
                      Point v2, 
                      Point v3);
    void AddPoint(vector<Point> point);
    //void AddPoint(double x, double y, double z);
    void GenerateTri(Point v1, Point v2, Point v3);
    void DeleteTri(int n, vector<Triangle> triangle);
    bool Iscirle(double x, double y, vector<Triangle> triangle);
    void DeleteFrame();


};

DelaunayTriangulation::DelaunayTriangulation(vector<Point> points){
    Point p_min = points[0];
    Point p_max = points[0];
    for (const auto & p :points){                       //遍历points上的节点
        if (p.x < p_min.x)  p_min.x = p.x;
        if (p.y < p_min.y)  p_min.y = p.y;
        if (p.z < p_min.z)  p_min.z = p.z;
        if (p.x > p_max.x)  p_max.x = p.x;
        if (p.y > p_max.y)  p_max.y = p.y;
        if (p.z > p_max.z)  p_max.z = p.z;
        //找出最边界的顶点设定为最小和最大的坐标点
        //构建包围盒，生成两个超三角形作为初始三角网络
        // GenerateTri(p1, p2, p3, p4);
    }
}
//析构函数实现
DelaunayTriangulation::~DelaunayTriangulation(){
    point.resize(0);
    triangle.resize(0);
}

//有参构造函数实现
void DelaunayTriangulation::Circumscribe(double &x_center, double &y_center, double &radius,  Point v1, Point v2,  Point v3){
    double x_1, x_2, x_3, y_1, y_2,y_3;
    x_1 = v1.x;
    y_1 = v1.y;
    x_2 = v2.x;
    y_2 = v2.y;
    x_3 = v3.x;
    y_3 = v3.y;
    x_center = (((y_2 - y_1)*(y_3*y_3 - y_1*y_1 + x_3*x_3 - x_1*x_1))
               -((y_3 - y_1)*(y_2*y_2 - y_1*y_1 + x_2*x_2 - x_1*x_1)))
               /((2*(x_3 - x_1)*(y_2 - y_1)) - (2*(x_2 - x_1)*(y_3 - y_1)));

    y_center = (((x_2 - x_1)*(x_3*x_3 - x_1*x_1 + y_3*y_3 - y_1*y_1))
               -((x_3 - x_1)*(x_2*x_2 - x_1*x_1 + y_2*y_2 - y_1*y_1)))
               /((2*(y_3 - y_1)*(x_2 - x_1)) - (2*(y_2 - y_1)*(x_3 - x_1)));

    radius = sqrt(pow(x_center - x_1,2) + pow(y_center - y_1, 2));

}
