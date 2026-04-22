
#include "include/mapbox/polylabel.hpp"

#include <fstream> 
#include <iostream>
#include <Windows.h>
#include <iomanip>
#include <stack>
#include "PolyCenter.hpp"
#include "io.h"
#include <sstream>
#include <string>
static int getrepeattimes() {
    int value;
    std::string input;
    while (true) {
        std::cout << "Please enter the number of repetitions £¨1-1000£©:";
        std::getline(std::cin, input);
        std::istringstream iss(input);
        if (iss >> value && iss.eof() && value >= 1 && value <= 1000) return value;
        return 10;
    }

}
static std::string GetExecutablePath() {
    char buffer[MAX_PATH];
    if (GetModuleFileNameA(NULL, buffer, MAX_PATH) > 0) {
        std::string path = buffer;
        size_t last_slash_idx = path.find_last_of("\\/");
        if (last_slash_idx != std::string::npos) {
            return path.substr(0, last_slash_idx);
        }
    }
    return "";
}
static std::string fixStringLength(const std::string& input, int length) {
    if (input.length() > length) {
        return input.substr(0, length);
    }
    else if (input.length() < length) {
        int padding = (length - input.length()) / 2;
        return std::string(padding, ' ') + input + std::string(length - input.length() - padding, ' ');
    }
    else {
        return input;
    }
}
static std::string fixDoubleLength(double value, int precision, int length, std::string suffix = "") {
    std::stringstream stream;
    stream << std::fixed << std::setprecision(precision) << value;
    return fixStringLength(stream.str() + suffix, length);
}
static std::vector< mapbox::geometry::polygon<double>> ReadforMapboxPolygons(const std::string& File) {
    std::vector< mapbox::geometry::polygon<double>> result;
    std::ifstream S(File, std::ios::binary);
    int numPolygons = 0;
    S.read(reinterpret_cast<char*>(&numPolygons), sizeof(numPolygons));
    for (int i = 0; i < numPolygons; ++i) {
        mapbox::geometry::linear_ring<double>  r;
        int numPoints = 0;
        S.read(reinterpret_cast<char*>(&numPoints), sizeof(numPoints));

        for (int j = 0; j < numPoints / 16; ++j) {
            double x = 0, y = 0;
            S.read(reinterpret_cast<char*>(&x), sizeof(x));
            S.read(reinterpret_cast<char*>(&y), sizeof(y));
            r.push_back(mapbox::geometry::point<double>(x, y));
        }
        mapbox::geometry::polygon<double> t = { r };
        result.push_back(t);
    }
    S.close();
    return result;
}
static mapbox::geometry::polygon<double> ReadforMapboxPolygon(const std::string& File) {
    mapbox::geometry::polygon<double> PG;
    std::ifstream S(File, std::ios::binary);
    int numLists = 0;
    S.read(reinterpret_cast<char*>(&numLists), sizeof(numLists));

    for (int i = 0; i < numLists; ++i) {
        mapbox::geometry::linear_ring<double>  r;
        int numPoints = 0;
        S.read(reinterpret_cast<char*>(&numPoints), sizeof(numPoints));

        for (int j = 0; j < numPoints / 16; ++j) {
            double x = 0, y = 0;
            S.read(reinterpret_cast<char*>(&x), sizeof(x));
            S.read(reinterpret_cast<char*>(&y), sizeof(y));
            r.push_back(mapbox::geometry::point<double>(x, y));
        }
        PG.push_back(r);
    }
    S.close();
    return PG;
}
static std::vector< std::vector<std::vector<PolyCenter::Vertex*>>> ReadforPolyCenterPolygons(const std::string& File) {
    std::vector < std::vector<std::vector<PolyCenter::Vertex*>>> R;
    std::ifstream S(File, std::ios::binary);
    if (!S) {
        std::cerr << "Cannot open file." << std::endl;
        return R;
    }

    int numLists = 0;
    S.read(reinterpret_cast<char*>(&numLists), sizeof(numLists));

    for (int i = 0; i < numLists; ++i) {
        std::vector<PolyCenter::Vertex*> A;
        int numPoints = 0;
        S.read(reinterpret_cast<char*>(&numPoints), sizeof(numPoints));

        for (int j = 0; j < numPoints / 16; ++j) {
            double x = 0, y = 0;
            S.read(reinterpret_cast<char*>(&x), sizeof(x));
            S.read(reinterpret_cast<char*>(&y), sizeof(y));
            A.push_back(new PolyCenter::Vertex(x, y));
        }
        std::vector<std::vector<PolyCenter::Vertex*>> B = { A };
        R.push_back(B);
    }
    S.close();
    return R;

}
static std::vector<std::vector<PolyCenter::Vertex*>> ReadforPolyCenter(const std::string& File) {
    std::vector<std::vector<PolyCenter::Vertex*>> R;
    std::ifstream S(File, std::ios::binary);

    if (!S) {
        std::cerr << "Cannot open file." << std::endl;
        return R;
    }

    int numLists = 0;
    S.read(reinterpret_cast<char*>(&numLists), sizeof(numLists));

    for (int i = 0; i < numLists; ++i) {
        std::vector<PolyCenter::Vertex*> A;
        int numPoints = 0;
        S.read(reinterpret_cast<char*>(&numPoints), sizeof(numPoints));
        for (int j = 0; j < numPoints / 16; ++j) {
            double x = 0, y = 0;
            S.read(reinterpret_cast<char*>(&x), sizeof(x));
            S.read(reinterpret_cast<char*>(&y), sizeof(y));
            A.push_back(new PolyCenter::Vertex(x, y));
        }
        R.push_back(A);
    }
    S.close();
    return R;

}
static std::vector<std::string> getAllFiles(std::string path, std::string fileType)
{
    long long hFile = 0;  struct _finddata_t fileinfo; std::string p;
    std::vector<std::string> r;
    if ((hFile = _findfirst(p.assign(path).append("\\*" + fileType).c_str(), &fileinfo)) != -1) {
        do {
            r.push_back(p.assign(path).append("\\").append(fileinfo.name));
        } while (_findnext(hFile, &fileinfo) == 0); 

        _findclose(hFile);
    }
    return r;
}

static void testpolylabel(mapbox::geometry::polygon<double> pts, int times, double precision, std::string s) {
    auto tick = GetTickCount64();  mapbox::geometry::point<double> c;
    for (int i = 0;i < times;i++) {
        c = mapbox::polylabel(pts, precision);
    }
    tick = GetTickCount64() - tick;
    auto radius = mapbox::detail::pointToPolygonDist(c, pts);
    std::cout << fixStringLength("Polylabel (" + s + ")", 24) << fixDoubleLength(c.x, 9, 18) << fixDoubleLength(c.y, 9, 18) << fixDoubleLength(radius, 9, 18) << fixDoubleLength(tick, 0, 16, "(ms)") << fixDoubleLength(static_cast<double>(tick) / times, 3, 16, "(ms)") << "\n";
}
static void testpolycircle(std::vector < std::vector<PolyCenter::Vertex*>> pts, int times) {
    auto tick = GetTickCount64(); auto c = PolyCenter::Circle();
    for (int i = 0;i < times;i++) {
        c = PolyCenter::Polygon(pts).maximumCircle();
    }
    tick = GetTickCount64() - tick;
    std::cout << fixStringLength("Polycenter", 24) << fixDoubleLength(c.centerx, 9, 18) << fixDoubleLength(c.centery, 9, 18) << fixDoubleLength(c.radius, 9, 18) << fixDoubleLength(tick, 0, 16, "(ms)") << fixDoubleLength(static_cast<double>(tick) / times, 3, 16, "(ms)") << "\n";
}

static std::string getname(std::string file) {
    std::string::size_type iPos = file.find_last_of('\\') + 1;
    std::string filename = file.substr(iPos, file.length() - iPos);
    return filename.substr(0, filename.rfind(".")) + ":";
}
static int count(std::vector < std::vector<PolyCenter::Vertex*>> pts) {
    int c = 0;
    for (auto& p : pts) {
        c += p.size();
    }
    return c;
}

int main () {
    std::cout << "This program is used to compare the performance of the Polycenter algorithm and the Polylabel algorithm. Please ensure that the test datas is in the same directory as this program. \n";
    std::cout << "\nIn order to accurately calculate the program operation time :\n";
    int times = getrepeattimes();
    std::cout << "\nThe program will repeat the operation " << times << " times.\n";

    std::cout << "\n" << fixStringLength("algorithm", 24) << fixStringLength("x", 18) << fixStringLength("y", 18) << fixStringLength("radius", 18) << fixStringLength("total time", 16) << fixStringLength("average time", 16) << "\n\n";
    auto path = GetExecutablePath();
    for (auto &ts : getAllFiles(path, "*.dat")) {
        auto pts1 = ReadforMapboxPolygon(ts);
        auto pts2 = ReadforPolyCenter(ts);
        std::cout << " " << getname(ts) << "(Number of vertices: " << count(pts2)<<  " )\n";
    
        testpolylabel(pts1, times, 1.0, "1.0");
        testpolylabel(pts1, times, 0.0001, "0.0001");
        testpolylabel(pts1, times, 0.00000001, "0.00000001");
        testpolycircle(pts2, times);
        std::cout << "\n";
    }
    std::cout << "Press any key to continue . . ." << std::endl;
    std::cin.ignore();
    std::cin.get();
    system("pause");
    return 0;
}

 
int main1 () {
    auto list1 = ReadforMapboxPolygons("D:\\USA.txt");
    auto list2 = ReadforPolyCenterPolygons("D:\\USA.txt");
    auto tick = GetTickCount64();
    tick = GetTickCount64();
    int i = 0;
    for (auto& l : list2) {
       PolyCenter::Polygon(l).maximumCircle();
    }

    std::cout << "PolyCenter:  " << GetTickCount64() - tick << "(ms)\n";

    for (auto& l : list1) {
        mapbox::polylabel(l, 1.0);
    }
    std::cout << "polylabel(1.0):  " << GetTickCount64() - tick << "(ms)\n";

    tick = GetTickCount64();
    for (auto& l : list1) {
        mapbox::polylabel(l, 0.0001);
    }
    std::cout << "polylabel(0.0001):  " << GetTickCount64() - tick << "(ms)\n";


    tick = GetTickCount64();
    for (auto& l : list1) {
        mapbox::polylabel(l, 0.00000001);
    }
    std::cout << "polylabel(0.00000001):  " << GetTickCount64() - tick << "(ms)\n";

 
    return 0;
}

 