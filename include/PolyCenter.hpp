#pragma once
#include <vector>
#include <deque>
#include <fstream>
#include <stack>
#include <algorithm>
namespace PolyCenter {
	constexpr double FloatACC = 0.0000000001, JudgeACC = 0.00000001;
	constexpr double PI = 3.14159265358979323846;
	static int CurrentArcNo = 0, CurrentRingNo = 0;
	static double maximumRadius = 0, maximumArea = 0, maximumCenterX = 0, maximumCenterY = 0;
	struct Circle {
		double centerx, centery, radius;   bool IsEmpty;
		Circle(double centerx, double centery, double radius) : centerx(centerx), centery(centery), radius(radius), IsEmpty(false) {}
		Circle() :IsEmpty(true), centerx(0), centery(0), radius(0) {}
	};
	struct Point2d {
		double x, y;
		Point2d(double x, double y) :x(x), y(y) {}
	};
	class Ring;
	class Vertex {
	public:
		Vertex(const double x, const double y) : x(x), y(y), m(0), n(0) {}
		Vertex* Above = nullptr, * Below = nullptr; Ring* ring = nullptr;
		double x, y, m, n;  bool concave = false;  //for point;
		double length = 0, sinL = 0, cosL = 0, sinA = 0, d = 0, cosA = 1;  //for Side
		double r1 = 0, r2 = 0, d1 = 0, d2 = 0; //for Bound
		int arc = 0;
		void setBelow(Vertex* next) {
			double dx = next->x - x;
			double dy = next->y - y;
			Below = next; next->Above = this;
			length = sqrt(dx * dx + dy * dy);
			sinL = dy / length;
			cosL = dx / length;
		}
		Point2d ArcCenter() const {
			double dt = maximumRadius - length / 2 < 0 ? 0 : sqrt(maximumRadius * maximumRadius - length * length / 4);
			return Point2d((x + Below->x) / 2 - dt / length * (Below->y - y), (y + Below->y) / 2 + dt / length * (Below->x - x));
		}
		void projectforLine(const double sinR, const  double cosR) {
			double m2 = Below->m, n2 = Below->n;
			cosA = cosL * cosR + sinL * sinR;
			sinA = sinL * cosR - cosL * sinR;
			if (fabs(n - n2) < 2 * JudgeACC) {
				cosA = 0; n = (n + n2) / 2;
				sinA = m < m2 ? 1 : -1;
			}
			else {
				if (fabs(m - m2) < 2 * JudgeACC) {
					sinA = 0; m = (m + m2) / 2;m2 = m;
					cosA = (n < n2) ? 1 : -1;
				}
				d = (m2 * n - m * n2) / (n - n2) * cosA;
				d1 = m * sinA / cosA + n;
				d2 = m2 * sinA / cosA + n2;
				r1 = fabs(m / cosA);
				r2 = fabs(m2 / cosA);
			}
		}
		void testRadiusforVertex(const double testL, double& testR) const {
			auto tp = sqrt((n - testL) * (n - testL) + m * m);
			if (tp < testR) testR = tp;
		}
		void testRadiusforSide(const double testL, double& testR) const {
			if ((testL - d2) * (d2 - d1) > 0 || (testL - d1) * (d1 - d2) > 0 || testR < fmin(r1, r2) + JudgeACC) return;
			auto tp = arc ? fabs(arc >= CurrentArcNo ? testR : maximumRadius - sqrt((testL - sinA) * (testL - sinA) + d * d)) : fabs(d + testL * sinA);
			if (tp < testR) testR = tp;
		}
		bool IsContactSide(const double testL, double testR) {//whether this side is in contact with the current circle.
			if (testR < fmin(r1, r2) - JudgeACC || (testL - d2) * (d2 - d1) > 0 || (testL - d1) * (d1 - d2) > 0) return false;
			return (testR - fabs(arc ? maximumRadius - sqrt((testL - sinA) * (testL - sinA) + d * d) : d + testL * sinA) > -JudgeACC);
		}
		void ChangeRing(Ring* newRing) {
			Vertex* current = this;
			do {
				current->ring = newRing;
				current = current->Below;
			} while (current != this);
		}
		Point2d Area() {
			Vertex* current = this; double areaMain = 0, areaArc = 0;
			do {
				areaMain += current->y * (current->Above->x - current->Below->x);
				if (current->arc) {
					double L = current->length;
					if (maximumRadius - L / 2 < FloatACC) {
						areaArc += maximumRadius * maximumRadius * PI / 2;
					}
					else {
						double dt = sqrt(maximumRadius * maximumRadius - L * L / 4);
						areaArc += maximumRadius * maximumRadius * asin(L / 2 / maximumRadius) - L * dt / 2;
					}
				}
				current->ring = ring;
				current = current->Below;
			} while (current != this);
			areaMain /= 2;
			if (abs(areaMain) < maximumRadius * JudgeACC) areaMain = 0;
			if (areaMain < 0) areaArc = -areaArc;
			return Point2d(areaMain, areaArc);
		}

		Point2d Middle(double sinR, double cosR) const {
			double x0 = 0, y0 = 0;
			if (arc) {
				double dt = maximumRadius - length / 2 < 0 ? 0 : sqrt(maximumRadius * maximumRadius - length * length / 4);
				x0 = (x + Below->x) / 2 + (dt - maximumRadius) * cosR;
				y0 = (y + Below->y) / 2 + (dt - maximumRadius) * sinR;
			}
			else {
				x0 = (x + Below->x) / 2; y0 = (y + Below->y) / 2;
			}
			return Point2d(x0, y0);
		}
	};

	namespace tools {
		std::vector<Vertex*> Vertexfordiscard;std::vector<Ring*> Ringfordiscard;
		static void clearrings() {
			for (auto& v : Vertexfordiscard) {
				std::vector<Vertex*> vs;
				auto f = v;
				do {
					vs.push_back(f); f = f->Below;
				} while (f != v);
				for (auto& t : vs) {
					delete(t);
				}
			}
			for (auto& r : Ringfordiscard) {
				delete(r);
			}
		}
		template<typename T1, typename T2>static double distance(const T1& A, const T2& B) {
			return sqrt((A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y));
		}
		template<typename T1, typename T2, typename T3>static double distance(const T1 V, const T2 V1, const  T3 V2) {
			return ((V.x - V1.x) * (V2.y - V1.y) - (V.y - V1.y) * (V2.x - V1.x)) / distance(V1, V2);
		}
		template <typename T> static void remove_element(std::vector<T>& vec, const T& value) {
			auto it = std::find(vec.begin(), vec.end(), value);
			if (it != vec.end()) {
				vec.erase(it);
			}
		}
		template<typename T1, typename T2, typename T3> static Point2d Per(const T1 p, const T2 P1, const T3 P2) {
			double rx = (P2.y - P1.y) * (P1.y - p.y) * (P2.x - P1.x) - (P2.y - P1.y) * (P2.y - P1.y) * P1.x - (P2.x - P1.x) * (P2.x - P1.x) * p.x;
			double ry = (P2.x - P1.x) * (P1.x - p.x) * (P2.y - P1.y) - (P2.x - P1.x) * (P2.x - P1.x) * P1.y - (P2.y - P1.y) * (P2.y - P1.y) * p.y;
			double rr = -(P2.y - P1.y) * (P2.y - P1.y) - (P2.x - P1.x) * (P2.x - P1.x);
			return Point2d(rx / rr, ry / rr);
		}
		template<typename T1, typename T2, typename T3> static bool atLeft(const T1 p, const T2 p1, const T3 p2, double err) {
			return ((p.x - p1.x) * (p2.y - p1.y) - (p.y - p1.y) * (p2.x - p1.x)) / distance(p1, p2) < err;
		}
		template<typename T1, typename T2, typename T3> static double Project(const T1 p, const T2 p1, const T3 p2) {
			return ((p.x - p1.x) * (p2.x - p1.x) + (p.y - p1.y) * (p2.y - p1.y)) / distance(p1, p2);
		}
		template<typename T1, typename T2> static Point2d GetCenter(T1* A, T2* B) {
			double l = distance(*A, *B);
			double dt = maximumRadius - l / 2 < 0 ? 0 : sqrt(maximumRadius * maximumRadius - l * l / 4);
			double x = (A->x + B->x) / 2 - (A->y - B->y) * dt / l;
			double y = (A->y + B->y) / 2 - (B->x - A->x) * dt / l;
			return Point2d(x, y);
		}
		static int PointPoint(double& testL, double& testR, Vertex* v1, Vertex* v2) {
			if (fabs(v2->n - v1->n) < JudgeACC) return 0;   //0 none 1 wrong -1 right 1
			double tl = (v2->m * v2->m + v2->n * v2->n - v1->m * v1->m - v1->n * v1->n) / (v2->n - v1->n) / 2;
			double tr = sqrt(v1->m * v1->m + (v1->n - tl) * (v1->n - tl));
			if ((tr - testR) * (testR + JudgeACC - tr) < 0) return -1;
			testR = tr;testL = tl;return 1;
		}
		static int PointLine(double& testL, double& testR, Vertex* V, Vertex* L) {
			double T1 = (L->Below->m - L->m) / L->length;
			double T2 = (L->m * (L->Below->n - L->n) - L->n * (L->Below->m - L->m)) / L->length;
			double A = T1 * T1 - 1;
			double B = 2 * T1 * T2 + 2 * V->n;
			double C = T2 * T2 - V->n * V->n - V->m * V->m;
			double tp = 0;
			if (fabs(A) < FloatACC) {
				tp = -C / B;
			}
			else {
				double det = sqrt(B * B - 4 * A * C);
				if (isnan(det) || det < 0) det = 0;
				double r1 = (-B + det) / 2 / A, r2 = (-B - det) / 2 / A;
				tp = (fabs(r1 - testL) < fabs(r2 - testL)) ? r1 : r2;
			}
			double tr = sqrt(V->m * V->m + (V->n - tp) * (V->n - tp));
			if ((tr - testR) * (testR + JudgeACC - tr) < 0) return -1;
			testR = tr; testL = tp;return 1;
		}
		static int LineLine(double& testL, double& testR, Vertex* L1, Vertex* L2) {
			double A = L2->length * (L1->Below->m - L1->m) - L1->length * (L2->Below->m - L2->m);
			double B = L1->length * (L2->m * (L2->Below->n - L2->n) - L2->n * (L2->Below->m - L2->m));
			double C = L2->length * (L1->m * (L1->Below->n - L1->n) - L1->n * (L1->Below->m - L1->m));
			if (fabs(A) < JudgeACC) return 0;
			double tl = (B - C) / A;
			double tr = fabs(((tl - L1->n) * (L1->Below->m - L1->m) + L1->m * (L1->Below->n - L1->n))) / L1->length;
			if ((tr - testR) * (testR + JudgeACC - tr) < 0) return -1;
			testR = tr;testL = tl;return 1;
		}
		static int PointCircle(double& testL, double& testR, Vertex* V, Vertex* Arc) {
			double rr = maximumRadius * maximumRadius;
			double dt = maximumRadius - Arc->length / 2 < 0 ? 0 : sqrt(rr - Arc->length * Arc->length / 4);
			double cm = (Arc->m + Arc->Below->m) / 2 + dt / Arc->length * (Arc->Below->n - Arc->n);
			double cn = (Arc->n + Arc->Below->n) / 2 - dt / Arc->length * (Arc->Below->m - Arc->m);
			double s = (rr + cn * cn + cm * cm - V->m * V->m - V->n * V->n) / 2;
			double A = rr - (V->n - cn) * (V->n - cn);
			double B = -2 * cn * rr - 2 * s * (V->n - cn);
			double C = rr * (cn * cn + cm * cm) - s * s;
			double tp = 0;
			if (fabs(A) < JudgeACC) {
				tp = -C / B;
			}
			else {
				double det = sqrt(B * B - 4 * A * C);
				if (isnan(det) || det < FloatACC) det = 0;
				double r1 = (-B + det) / 2 / A, r2 = (-B - det) / 2 / A;
				tp = ((r1 - Arc->Below->n) * (cm - Arc->Below->m) + Arc->Below->m * (cn - Arc->Below->n) < 0 &&
					(r1 - Arc->n) * (cm - Arc->m) + Arc->m * (cn - Arc->n) > 0) ? r1 : r2;
			}
			double tr = sqrt(V->m * V->m + (V->n - tp) * (V->n - tp));
			if ((tr - testR) * (testR + JudgeACC - tr) < 0) return -1;
			testR = tr; testL = tp;return 1;
		}
		static int LineCircle(double& testL, double& testR, Vertex* L, Vertex* Arc) {
			double dt = maximumRadius - Arc->length / 2 < 0 ? 0 : sqrt(maximumRadius * maximumRadius - Arc->length * Arc->length / 4);
			double cm = (Arc->m + Arc->Below->m) / 2 + dt / Arc->length * (Arc->Below->n - Arc->n);
			double cn = (Arc->n + Arc->Below->n) / 2 - dt / Arc->length * (Arc->Below->m - Arc->m);
			double T = (L->Below->m - L->m) / L->length;
			double S = maximumRadius + (L->m * (L->Below->n - L->n) - L->n * (L->Below->m - L->m)) / L->length;
			double A = T * T - 1;
			double B = 2 * (T * S + cn);
			double C = S * S - cn * cn - cm * cm;
			double tp = 0;
			if (fabs(A) < JudgeACC) {
				tp = -C / B;
			}
			else {
				double det = sqrt(B * B - 4 * A * C);
				if (isnan(det) || det < JudgeACC) det = 0;
				double r1 = (-B + det) / 2 / A, r2 = (-B - det) / 2 / A;
				tp = (r1 - L->n) * (L->m - L->Below->m) + L->m * (L->n - L->Below->n) > 0 ? r1 : r2;
			}
			double tr = fabs((tp - L->n) * (L->m - L->Below->m) + L->m * (L->n - L->Below->n)) / L->length;
			if ((tr - testR) * (testR + JudgeACC - tr) < 0) return -1;
			testR = tr;testL = tp;return 1;
		}
		static int CircleCircle(double& testL, double& testR, Vertex* Arc1, Vertex* Arc2) {
			double dt1 = maximumRadius - Arc1->length / 2 < 0 ? 0 : sqrt(maximumRadius * maximumRadius - Arc1->length * Arc1->length / 4);
			double cm1 = (Arc1->m + Arc1->Below->m) / 2 + dt1 / Arc1->length * (Arc1->Below->n - Arc1->n);
			double cn1 = (Arc1->n + Arc1->Below->n) / 2 - dt1 / Arc1->length * (Arc1->Below->m - Arc1->m);
			double dt2 = maximumRadius - Arc2->length / 2 < 0 ? 0 : sqrt(maximumRadius * maximumRadius - Arc2->length * Arc2->length / 4);
			double cm2 = (Arc2->m + Arc2->Below->m) / 2 + dt2 / Arc2->length * (Arc2->Below->n - Arc2->n);
			double cn2 = (Arc2->n + Arc2->Below->n) / 2 - dt2 / Arc2->length * (Arc2->Below->m - Arc2->m);
			double tl = (cm2 * cm2 + cn2 * cn2 - cm1 * cm1 - cn1 * cn1) / (cn2 - cn1) / 2;
			double tr = maximumRadius - sqrt(cm1 * cm1 + (cn1 - tl) * (cn1 - tl));
			if ((tr - testR) * (testR + JudgeACC - tr) < 0) return -1;
			testR = tr;testL = tl;return 1;
		}


	}
	using namespace tools;

	class Position {
	public:
		Vertex* base; double x, y; double note = 0;
		Position(Vertex* v) : base(v), x(v->x), y(v->y) {}
		Position(Vertex* v, double x, double y) : base(v), x(x), y(y) {
			if (distance(*this, *v->Below) < JudgeACC)  base = v->Below;
		}
		Vertex*& InsertVertex() { //Create a new vertex at the current position and insert it into the queue.
			if (distance(*this, *base) < JudgeACC) return base;
			auto v = new Vertex(x, y);
			v->ring = base->ring;v->arc = base->arc;
			auto bl = base->Below;
			base->setBelow(v); v->setBelow(bl);
			base = v;v->concave = true;
			return base;
		}
	};
	class Ring {
	public:
		double allowRadius = 0; bool notrim = false; bool isOut = false; bool tag = false;
		Point2d area = Point2d(0, 0); //Main area(x) and arc area(y); a mainarea < 0 indicates a hole.
		Vertex* first = nullptr;
		Ring(std::vector <Vertex*>& Vs) {
			double areaMain = 0, MAX = 0;
			first = Vs.front();auto Last = first;
			for (size_t i = 1; i < Vs.size(); ++i) {
				if (distance(*Vs[i], *Last) < JudgeACC) continue;
				if (abs(Vs[i]->x) > MAX) MAX = abs(Vs[i]->x);
				if (abs(Vs[i]->y) > MAX) MAX = abs(Vs[i]->y);
				Last->setBelow(Vs[i]);
				Last = Vs[i];
			}
			//ACCURACY = MAX * 0.000000000001;
			if (distance(*Last, *first) < JudgeACC) { first = first->Below; }
			Last->setBelow(first);
			Last = first;
			do {
				areaMain += Last->y * (Last->Above->x - Last->Below->x);
				concavity(Last);
				Last->arc = 0; Last->ring = this;
				Last = Last->Below;
			} while (Last != first);
			areaMain /= 2;
			area = Point2d(areaMain, 0);
			isOut = area.x > 0;
		}
		inline static void concavity(Vertex* v) {
			v->concave = ((v->Below->x - v->Above->x) * (v->y - v->Above->y) - (v->Below->y - v->Above->y) * (v->x - v->Above->x)) /
				(v->Above->length * v->Above->length) > -JudgeACC;
		}
		Ring(Vertex* first) :first(first) { first->ring = this; area = first->Area(); }
		bool Contains(Vertex* Pt) const {// whether a point is inside the ring.
			bool inside = false; Vertex* F = first;
			do {
				if ((F->y > Pt->y) != (F->Below->y > Pt->y)) {
					if (F->arc && distance(*Pt, F->ArcCenter()) < maximumRadius && atLeft(*Pt, *F->Below, *F, JudgeACC)) {
						return true;
					}
					double slope = (F->Below->x - F->x) * (Pt->y - F->y) / (F->Below->y - F->y);
					if (Pt->x < slope + F->x) inside = !inside;
				}
				F = F->Below;
			} while (F != first);
			return inside;
		}

		std::vector<Ring*> trim(std::vector<Ring*>& Holes) {
			std::vector<part> r;std::vector<Ring*> result;
			if (allowRadius >= maximumRadius - JudgeACC) { result.push_back(this);return result; }
			arctrim(r, first->Below, true, first->Below, first);
			for (auto it = Holes.begin();it < Holes.end();) {
				if (arctrim(r, (*it)->first, false, first->Below, first)) it = Holes.erase(it);
				else it++;
			}
			if (r.size() == 1) {
				auto first = splitout(r[0].fromV, r[0].toV, true);
				if (first) {
					bool isout = first->ring->isOut;
					auto w = new Ring(first); w->isOut = isout;
					result.push_back(w);
				}
				return result;
			}
			std::sort(r.begin(), r.end(), [](const part& A, const part& B) {return A.note1 < B.note1;});
			std::vector<part> out, in;
			for (auto& t : r) {
				if (t.note1 < t.note2) out.push_back(t);
				else(in.push_back(t));
			}
			std::sort(out.begin(), out.end(), [](const part& A, const part& B) {return A.note1 < B.note1;});
			std::sort(in.begin(), in.end(), [](const part& A, const part& B) {return A.note2 < B.note2;});
			int currentin = 0;
			for (auto& t : out) {
				auto ft = splitout(t.fromV, t.toV, true);
				if (!ft) continue;
				auto r = new Ring(ft);
				Vertex* v1 = r->first->Below, * v2 = r->first;
				if (r->area.x + r->area.y > maximumArea) {
					bool flag = false;
					for (auto it = in.begin();it < in.end();) {
						if (it->note2 > t.note1) {
							if (it->note1 > t.note2) break;
							auto ff = splitout(it->fromV, it->toV, true);
							if (!ff) break;
							auto l = Ring(ff);
							Vertex* v3 = l.first->Below, * v4 = l.first; double L = v4->length;
							close(v1, v4);close(v3, v2);
							t.note1 = it->note1; v1 = v3; flag = true;
							it = in.erase(it);
						}
						else {
							it++;
						}
					}
					if (flag) r->area = r->first->Area();
					result.push_back(r);
				}
				else {
					Vertexfordiscard.push_back(r->first);
				}
			}
			for (auto it = in.begin();it < in.end(); it++) {
				Vertexfordiscard.push_back(it->fromV);
			}
			return result;
		}
		static void close(Vertex* f, Vertex* t) {
			t->setBelow(f); t->arc = CurrentArcNo;
			f->concave = true;t->concave = true;
		}
		static Vertex* splitout(Vertex* fromV, Vertex* toV, bool discardrest) {
			if (fromV->Below == toV) {
				if (discardrest) Vertexfordiscard.push_back(fromV);
				return nullptr;
			}
			if (discardrest && fromV->ring == toV->ring) {
				if (fromV->Above != toV) {
					auto f = fromV->Above; f->setBelow(toV->Below);
					Vertexfordiscard.push_back(f);
					close(fromV, toV);
				}
				return toV;
			}
			else {
				if (fromV->ring != toV->ring) {
					if (fromV->ring->isOut) {
						toV->ChangeRing(fromV->ring);
					}
					else {
						fromV->ChangeRing(toV->ring);
					}
				}
				auto f = new Vertex(fromV->x, fromV->y);
				auto t = new Vertex(toV->x, toV->y);
				f->setBelow(fromV->Below);
				f->arc = fromV->arc;
				toV->Above->setBelow(t);
				f->ring = fromV->ring;
				t->ring = fromV->ring;
				close(f, t);
				close(toV, fromV);
				return t;
			}
		}

		int getClipper() {
			if (continuousarc()) { area = first->Area(); }
			Vertex* f = first->Below, * t = first;
			if (f->Below == t) return 0;
			if (approximateParallel(f, t)) return 1;

			double acc = fmax(JudgeACC, sqrt(2 * maximumRadius * JudgeACC));
			while (true) {
				center = GetCenter(f, t);
				auto ff = ontooutside(f, t, acc);
				if (!ff) return 0;
				auto tt = ontooutsiderev(f, t, acc);
				if (!tt) return 0;
				if (ff == f && tt == t) break;
				f = ff;t = tt;
				if (atLeft(center, *t, *f, acc)) break;
			}
			if (f != first->Below) f = ontooutrev(first->Below, f);
			if (t != first) t = ontoout(t, first);
			if (f == nullptr) f = first->Below;
			if (t == nullptr) t = first;
			if (f == first->Below && t == first) return 1;
			first = splitout(f, t, true);
			if (!first) return 0;
			area = first->Area(); return 2;
		}
	private:
		bool continuousarc() {
			bool flag = false;
			if (first->Below->length < 2 * JudgeACC) removeVertex(first->Below);
			if (first->Above->length < 2 * JudgeACC) first = removeVertex(first);
			if (first->arc) {
				if (first->Above->arc && radius3point(first->Below, first, first->Above) < maximumRadius + JudgeACC) {
					first = removeVertex(first);  allowRadius = 0;flag = true;
				}
				else if (first->Below->arc && radius3point(first, first->Below, first->Below->Below) < maximumRadius + JudgeACC) {
					removeVertex(first->Below);allowRadius = 0;flag = true;
				}
			}
			if (flag) area = first->Area();
			return flag;
		}
		double radius3point(Vertex* v1, Vertex* v2, Vertex* v3) {
			double  a = sqrt((v2->x - v3->x) * (v2->x - v3->x) + (v2->y - v3->y) * (v2->y - v3->y));
			double  b = sqrt((v1->x - v3->x) * (v1->x - v3->x) + (v1->y - v3->y) * (v1->y - v3->y));
			double  c = sqrt((v1->x - v2->x) * (v1->x - v2->x) + (v1->y - v2->y) * (v1->y - v2->y));
			double  d = v1->x * (v2->y - v3->y) + v2->x * (v3->y - v1->y) + v3->x * (v1->y - v2->y);
			return a * b * c / 2 / abs(d);
		}
		static Vertex* removeVertex(Vertex* v) {
			auto t = v->Above;
			t->setBelow(v->Below);
			v->setBelow(v); Vertexfordiscard.push_back(v);
			return t;
		}
		bool approximateParallel(Vertex*& f, Vertex*& t) {

			double acc = 2 * sqrt(2 * maximumRadius * JudgeACC);
			if (fabs(distance(*t, *f) - 2 * maximumRadius) > 2 * JudgeACC) return false;
			double x0 = (f->x + t->x) / 2, y0 = (f->y + t->y) / 2;
			double sinR = (f->x - t->x) / t->length, cosR = (t->y - f->y) / t->length;
			do {
				f = f->Below;
				f->m = (f->y - y0) * cosR - (f->x - x0) * sinR;
			} while (fabs(f->m + maximumRadius) < f->Above->length / maximumRadius * acc);
			do {
				t = t->Above;
				t->m = (t->y - y0) * cosR - (t->x - x0) * sinR;
			} while (fabs(t->m - maximumRadius) < t->length / maximumRadius * acc);

			f = f->Above;t = t->Below;
			f->n = (f->x - x0) * cosR + (f->y - y0) * sinR;
			t->n = (t->x - x0) * cosR + (t->y - y0) * sinR;

			if (f->Below == t || f->Below == t->Above) return false;
			if (fabs(f->n - t->n) < 4 * JudgeACC) {

			}
			else if (f->n > t->n) {
				auto l = Project(*t, *f->Above, *f);
				if (f->length > acc) {
					f = f->Above;
					f->x += f->cosL * (l + acc);
					f->y += f->sinL * (l + acc);
					f->length -= l + acc;
				}
			}
			else {
				auto l = Project(*f, *t->Below, *t);
				if (t->length > acc) {
					t = t->Below;
					t->x -= t->Above->cosL * (l + acc);
					t->y -= t->Above->sinL * (l + acc);
					t->Above->length -= l + acc;
				}
			}
			first = splitout(f, t, true);
			return true;
		}

		struct part
		{
			Vertex* fromV, * toV; double note1, note2;
			part(Vertex* toV, Vertex* fromV, Vertex* arcfrom) : fromV(fromV), toV(toV), note1(distance(*arcfrom, *fromV)), note2(distance(*arcfrom, *toV)) {}

		};
		Point2d center = Point2d(0, 0);//, arcfrom = Point2d(0, 0), arcto = Point2d(0, 0);
		bool arctrim(std::vector<part>& result, Vertex* f, bool self, Vertex* arcfrom, Vertex* arcto) const {

			int current = -2;
			Point2d p1(0, 0), p2(0, 0);
			if (!self) {
				do {
					f = f->Below; double temp = distance(*f, center) - maximumRadius;
					current = temp > JudgeACC ? 1 : temp < -JudgeACC ? -1 : current < 0 ? -2 : 0;
				} while (current == 0);
			}
			auto t = self ? f->Above : f;
			std::deque<Vertex*>list;
			do {
				auto tp = Intersection(true, current, f, p1, p2);
				if (tp) {
					if (atLeft(p1, *arcfrom, *arcto, JudgeACC)) {
						f = Position(f, p1.x, p1.y).InsertVertex();
						list.push_back(f);
					}
					if (tp == 2 && atLeft(p2, *arcfrom, *arcto, JudgeACC)) {
						f = Position(f, p2.x, p2.y).InsertVertex();
						list.push_back(f);
					}
				}
				f = f->Below;
			} while (f != t);
			if (self) {
				if (list.empty() || list.front() != t->Below) list.push_front(t->Below);
				if (list.back() != t) {
					if ((list.size() & 1) == 0) {
						list.back() = t;
					}
					else {
						list.push_back(t);
					}
				}
			}
			else {
				if (list.size() == 0) return false;
				f = list[0]->Below;
				if (distance(*f, center) - maximumRadius > JudgeACC && atLeft(*f, *arcfrom, *arcto, JudgeACC))  std::rotate(list.begin(), list.end() - 1, list.end());
			}
			for (auto it = list.begin();it < list.end(); ) {
				auto a = *it++, b = *it++;
				result.emplace_back(b, a, arcfrom);
			}
			return true;
		}
		Vertex* ontooutsiderev(Vertex* f, Vertex* t, double acc) const {
			int current = -2; Point2d p1(0, 0), p2(0, 0); auto tt = t;
			while (t != f) {
				int w = Intersection(false, current, t, p1, p2);
				if (w) {
					if (atLeft(p1, *tt, *f, acc)) return Position(t->Above, p1.x, p1.y).InsertVertex();
					if (w == 2 && atLeft(p1, *tt, *f, acc)) return Position(t->Above, p2.x, p2.y).InsertVertex();
				}
				t = t->Above;
			}
			return nullptr;
		}
		Vertex* ontooutside(Vertex* f, Vertex* t, double acc) const {
			int current = -2; Point2d p1(0, 0), p2(0, 0); auto ff = f;
			while (f != t) {
				int w = Intersection(true, current, f, p1, p2);
				if (w) {
					if (atLeft(p1, *t, *ff, acc)) return Position(f, p1.x, p1.y).InsertVertex();
					if (w == 2 && atLeft(p1, *t, *ff, acc)) return Position(f, p2.x, p2.y).InsertVertex();
				}
				f = f->Below;
			}
			return nullptr;
		}
		Vertex* ontooutrev(Vertex* f, Vertex* t) const {
			int current = -2; Point2d p1(0, 0), p2(0, 0); auto tt = t;
			while (t != f) {
				if (Intersection(false, current, t, p1, p2)) return Position(t->Above, p1.x, p1.y).InsertVertex();
				t = t->Above;
			}
			return nullptr;
		}
		Vertex* ontoout(Vertex* f, Vertex* t) const {
			int current = -2; Point2d p1(0, 0), p2(0, 0); auto ff = f;
			while (f != t) {
				if (Intersection(true, current, f, p1, p2)) return Position(f, p1.x, p1.y).InsertVertex();
				f = f->Below;
			}
			return nullptr;
		}

		int Intersection(bool anti, int& current, Vertex* v1, Point2d& p1, Point2d& p2) const { //inter count
			Vertex* v2 = anti ? v1->Below : v1->Above;Vertex* v = anti ? v1 : v2;
			double length = v->length; double temp = distance(*v2, center) - maximumRadius;
			int next = temp > JudgeACC ? 1 : temp < -JudgeACC ? -1 : current < 0 ? -2 : 0;
			if (next == 0) {// out - on
				double L2 = ((center.x - v1->x) * (v2->x - v1->x) + (center.y - v1->y) * (v2->y - v1->y)) / v->length;
				if (L2 < v->length - JudgeACC) {
					double L1 = fabs((center.x - v1->x) * (v2->y - v1->y) - (center.y - v1->y) * (v2->x - v1->x)) / v->length;
					double h = sqrt(maximumRadius * maximumRadius - L1 * L1);current = -2;
					p1 = Point2d(v1->x + (L2 - h) * (v2->x - v1->x) / v->length, v1->y + (L2 - h) * (v2->y - v1->y) / v->length);
					return 1;
				}
				else {
					current = 2; return 0;
				}
			}
			int cuut = current;current = next;
			if (cuut * next == -1) {
				if (v->arc) {
					auto c = v->ArcCenter(); double d = distance(c, center);
					double h = sqrt(maximumRadius * maximumRadius - d * d / 4);
					p1 = Point2d((center.x + c.x) / 2 - h * (center.y - c.y) / d, (center.y + c.y) / 2 + h * (center.x - c.x) / d);
					if (atLeft(p1, *v, *v->Below, JudgeACC))
						p1 = Point2d((center.x + c.x) / 2 + h * (center.y - c.y) / d, (center.y + c.y) / 2 - h * (center.x - c.x) / d);
				}
				else {
					double L1 = ((center.x - v1->x) * (v2->y - v1->y) - (center.y - v1->y) * (v2->x - v1->x)) / v->length;
					double L2 = ((center.x - v1->x) * (v2->x - v1->x) + (center.y - v1->y) * (v2->y - v1->y)) / v->length;
					double h = next * sqrt(maximumRadius * maximumRadius - L1 * L1);
					p1 = Point2d(v1->x + (L2 + h) * (v2->x - v1->x) / v->length, v1->y + (L2 + h) * (v2->y - v1->y) / v->length);
				}
				return 1;
			}
			if (cuut == 2 && next == -1) {
				p1 = Point2d(v1->x, v1->y); return 1;
			}
			if (next == 1 && v->arc == 0) {
				double L1 = fabs((center.x - v1->x) * (v2->y - v1->y) - (center.y - v1->y) * (v2->x - v1->x)) / v->length;
				if (cuut == -2) {
					double L2 = ((center.x - v1->x) * (v2->x - v1->x) + (center.y - v1->y) * (v2->y - v1->y)) / v->length;
					if (L2 < JudgeACC) {
						p1 = Point2d(v1->x, v1->y);
					}
					else {
						double h = maximumRadius < L1 ? 0 : sqrt(maximumRadius * maximumRadius - L1 * L1);
						p1 = Point2d(v1->x + (L2 + h) * (v2->x - v1->x) / v->length, v1->y + (L2 + h) * (v2->y - v1->y) / v->length);
					}
					return 1;
				}
				else if (cuut == 1) {
					if (L1 < maximumRadius + JudgeACC) {
						double L2 = ((center.x - v1->x) * (v2->x - v1->x) + (center.y - v1->y) * (v2->y - v1->y)) / v->length;
						if (L2 > JudgeACC && L2 < v->length - JudgeACC) {
							double h = sqrt(maximumRadius * maximumRadius - L1 * L1);
							p1 = Point2d(v1->x + (L2 - h) * (v2->x - v1->x) / v->length, v1->y + (L2 - h) * (v2->y - v1->y) / v->length);
							p2 = Point2d(v1->x + (L2 + h) * (v2->x - v1->x) / v->length, v1->y + (L2 + h) * (v2->y - v1->y) / v->length);
							return 2;
						}
					}
				}
			}
			if (next == 1 && v->arc && cuut == -2) { p1 = Point2d(v1->x, v1->y); return 1; }
			return 0;
		}
	};
	class IntervalsforL {
	public:
		IntervalsforL(std::vector<double> s) : s(s) {}
		double empty() {
			return s.empty();
		}
		bool isNarrow() {
			return s.empty() || s.back() - s.front() < JudgeACC;
		}
		double Middle() {
			double last = -1.0, r = 0.0;
			for (auto it = s.begin();it != s.end();) {
				double a = *it++, b = *it++;
				if (b - a > last) {
					last = b - a; r = (a + b) / 2;
				}
			}
			return r;
		}
		Point2d Max() {
			double last = -1.0; Point2d pt(0, 0);
			for (auto it = s.begin();it != s.end();) {
				double a = *it++, b = *it++;
				if (b - a > last) {
					last = b - a; pt = Point2d(a, b);
				}
			}
			return pt;
		}
		bool isVertexUseful(Vertex* v, double& Radius) const {
			auto x = v->n - Radius, y = v->n + Radius, mm = v->m * v->m;
			for (auto it = s.begin(); it != s.end(); ) {
				if ((*it++ - y) * (*it++ - x) + mm < 4 * JudgeACC) return true;
			}
			return false;
		}
		bool isSideUseful(Vertex* v, double& Radius) const {
			int i = Pos(v->d1);
			return ((i & 1) == 1 || i != Pos(v->d2)) && fmin(v->r1, v->r2) < (s.back() - s.front()) / 2.0 + Radius + 2 * JudgeACC;
		}
		void vertexConstraint(const Vertex* v, const double R) {
			if (R - fabs(v->m) > JudgeACC) {
				double det = sqrt(R * R - v->m * v->m);
				Minus(v->n - det, v->n + det);
			}

		}
		void sideConstraint(const Vertex* v, const double R) {
			if (v->arc) { arcConstraint(v, R); return; }
			if (R < fmin(v->r1, v->r2) + JudgeACC) return;
			auto down = fmin(v->d1, v->d2), up = fmax(v->d1, v->d2);
			if (R > fmax(v->r1, v->r2) + JudgeACC) {
				Minus(down, up);return;
			}
			double temp = -(v->d + R) / v->sinA;
			if (v->sinA < 0) {
				if (v->r2 > v->r1) {
					Minus(v->d1, temp);
				}
				else {
					Minus(v->d2, temp);
				}
			}
			else {
				if (v->r2 > v->r1) {
					Minus(temp, v->d1);
				}
				else {
					Minus(temp, v->d2);
				}
			}
		}
		void arcConstraint(const Vertex* v, const double R) {
			double thisR = (v->arc >= CurrentArcNo) ? R : maximumRadius; //sub-polygon's original arcs: maximumRadius; process arcs: R. same below.
			if (R < fmin(v->r1, v->r2) + JudgeACC)  return;
			auto down = fmin(v->d1, v->d2), up = fmax(v->d1, v->d2);
			double cm = v->d, cn = v->sinA;
			if ((down - cn) * (cn - up) > 0) {
				if (R > thisR - fabs(cm) + JudgeACC) {
					Minus(down, up);
				}
				else if (R > fmax(v->r1, v->r2) - JudgeACC) {
					double temp = thisR - R - fabs(cm) < 0 ? 0 : sqrt((thisR - R) * (thisR - R) - cm * cm);
					if (down < cn - temp)  Minus(down, cn - temp);
					if (cn + temp < up) Minus(cn + temp, up);
				}
				else if (R > fmin(v->r1, v->r2) - JudgeACC) {
					double temp = thisR - R - fabs(cm) < 0 ? 0 : sqrt((thisR - R) * (thisR - R) - cm * cm);
					if ((v->d1 - v->d2) * (v->r1 - v->r2) > 0) {
						Minus(down, cn - temp);
					}
					else {
						Minus(cn + temp, up);
					}
				}
			}
			else {
				if (R > fmax(v->r1, v->r2) - JudgeACC) {
					Minus(down, up);
				}
				else {
					double temp = thisR - R - fabs(cm) < 0 ? 0 : sqrt((thisR - R) * (thisR - R) - cm * cm);
					if ((v->d1 - v->d2) * (v->r1 - v->r2) > 0) {
						Minus(down - JudgeACC, cn - temp);
					}
					else {
						Minus(cn + temp, up + JudgeACC);
					}
				}

			}
		}
		void Minus(double lower, double upper) {
			int i = -1, j = s.size();
			for (auto a = 0; a < s.size(); ++a) {
				if (lower < s[a]) {
					i = a;break;
				}
			}
			if (i == -1) return;
			for (int b = s.size() - 1; b >= 0; --b) {
				if (upper > s[b]) {
					j = b;
					break;
				}
			}
			if (j == s.size()) return;
			s.erase(s.begin() + i, s.begin() + j + 1);
			if ((j & 1) == 0) s.insert(s.begin() + i, upper);
			if ((i & 1) == 1) s.insert(s.begin() + i, lower);

		}
		bool Contains(double d) {
			return Pos(d) & 1;
		}
	private:
		std::vector<double> s;
		int Pos(const double current) const {
			for (auto a = 0; a < s.size(); ++a) {
				if (current < s[a]) {
					return   a;
				}
			}
			return s.size();
		}

	};

	class Datumline {
	public:
		double x0, y0, sinR, cosR, radius = 0;
		Datumline(std::deque<Ring*>& rings) :rings(rings) {
			auto f = rings.front()->first;
			double x1 = f->x, y1 = f->y, x2 = f->Below->x, y2 = f->Below->y, l = f->length;
			cosR = (y1 - y2) / l;sinR = (x2 - x1) / l;
			auto pt = f->Middle(sinR, cosR);
			x0 = pt.x;y0 = pt.y;
		}
		std::vector<Ring*> GetPiece() {
			std::vector<Ring*> result; bool IsBound = false; double area = rings.front()->area.x;
			do {
				double testL = 0, lasttestL = -100;
				auto pp = InitFilter();
				while (fabs(lasttestL - testL) > FloatACC && !pp.isNarrow()) {
					if (radius > 0) Constraint(pp, radius);
					lasttestL = testL; testL = pp.Middle();
					radius = GetRadius(testL);
				}
				IsBound = split(pp, result);
				if (rings.size() && abs(rings.front()->area.x - area) < JudgeACC) { //rings.size() && abs(rings.front()->area.x - area) < JudgeACC
					auto r = rings.front();r->notrim = true; r->first = r->first->Below;
				}
			} while (!IsBound);

			return result;
		}
	private:

		bool split(IntervalsforL& pp, std::vector<Ring*>& result) {
			double firstarea = rings.front()->area.x;
			double testL = 0; bool isNormal = false;
			auto p = getPreciseStadiumContacts(pp, testL, isNormal);
			x0 += testL * cosR;y0 += testL * sinR;
			if (radius > maximumRadius) {
				maximumRadius = radius; maximumArea = PI * (radius + JudgeACC) * (radius + JudgeACC);
				maximumCenterX = x0;maximumCenterY = y0;
			}
			Point2d center(x0, y0);
			if (p.size() == 0) {
				rings.front()->first = rings.front()->first->Below; return true;
			}
			else if (p.size() == 1) {
				double t = sinR; sinR = cosR;cosR = -t;
				auto main = p.front().base;
				main->ring->first = main;
				return false;
			}
			else {
				Vertex* mainfirst = nullptr;
				for (auto it = p.begin();it != p.end();) {
					removering(it++->base->ring);
				}
				bool haschanged = false;
				for (auto it = p.rbegin(); it != p.rend();) {
					auto t = it++, f = (it == p.rend()) ? p.rbegin() : it;
					int s = isNormal ? 1 : sideflag(*f, *t);
					if (s == 0) {
						if (!haschanged) {
							mainfirst = t->InsertVertex(); haschanged = true;
							double tt = sinR; sinR = cosR;cosR = -tt;
						}
					}
					else if (s == 1) {
						bool fg = forNext(*f, *t);
						auto ff = f->InsertVertex(), tt = t->InsertVertex();
						bool isOut = ff->ring->isOut || tt->ring->isOut;
						bool isJoin = ff->ring != tt->ring;
						auto first = Ring::splitout(ff, tt, false);
						if (!first) continue;
						if (fg) mainfirst = first;
						if (!isJoin) {
							auto nR = new Ring(first); nR->isOut = isOut;

							result.push_back(nR);
						}
					}
				}
				if (mainfirst) {
					bool flag = false;
					if (mainfirst->arc && (mainfirst->Below->arc || mainfirst->Above->arc)) flag = true;
					remove_element(result, mainfirst->ring);
					if (flag) mainfirst->ring->area = mainfirst->Area();
					mainfirst->ring->first = mainfirst;
					rings.push_back(mainfirst->ring);
					return flag;
				}
				else {
					return true;
				}

			}



		}
		static int sideflag(Position f, Position t) {
			int i = sameside(f, t);
			return i == -1 ? 0 : (i == 1 || i == 3 ? 1 : 2);
		}
		static int sameside(Position& f, Position& t) {
			for (int i = 0; i < 4;i++) {
				if (isin(f, i) && isin(t, i)) return i;
			}
			return -1;
		}
		static bool isin(Position& v, int i) {
			if (i == 3 && v.note < JudgeACC) return true;
			return v.note > i - JudgeACC && v.note < i + 1 + JudgeACC;
		}

		bool forNext(Position f, Position t) {
			if (atLeft(Point2d(x0, y0), f, t, fmax(JudgeACC, sqrt(2 * maximumRadius * JudgeACC)))) return false;
			int tag = (distance(*f.base, f) < JudgeACC ? 0 : f.base->arc ? 20 : 10) + (distance(*t.base, t) < JudgeACC ? 0 : t.base->arc ? 2 : 1);
			bool flag = false;Point2d pt(0, 0);
			switch (tag)
			{
			case 2:
				flag = DirectionforArcPoint(pt, t.base, f.base); break;
			case 12:
				flag = DirectionforArcLine(pt, f.base, t.base); break;
			case 20:
				flag = DirectionforArcPoint(pt, f.base, t.base); break;
			case 21:
				flag = DirectionforArcLine(pt, f.base, t.base); break;
			case 22:
				flag = DirectionforArcArc(pt, f.base, t.base);break;
			}
			if (flag) {
				double l = distance(pt, Point2d(x0, y0));
				sinR = (pt.y - y0) / l;cosR = (pt.x - x0) / l;
			}
			else {
				double l = distance(f, t);
				sinR = (f.x - t.x) / l;cosR = (t.y - f.y) / l;
			}
			return true;
		}
		bool DirectionforArcPoint(Point2d& dir, const Vertex* arc, const  Vertex* pt) {
			auto center = arc->ArcCenter(); double l = distance(*pt, center);
			if (abs(l) < 4 * JudgeACC) return false;
			l = (maximumRadius + l) / 2 / l;
			dir = Point2d(pt->x + (center.x - pt->x) * l, pt->y + (center.y - pt->y) * l);
			return true;
		}
		bool DirectionforArcLine(Point2d& dir, const Vertex* arc, const Vertex* line) {
			auto center = arc->ArcCenter(); double l = distance(center, *line, *line->Below);
			if (abs(l) < 4 * JudgeACC) return false;
			l = (maximumRadius + l) / 2;
			dir = Point2d(center.x - l * line->sinL, center.y + l * line->cosL);
			return true;
		}
		bool DirectionforArcArc(Point2d& dir, const Vertex* arc1, const Vertex* arc2) {
			auto center1 = arc1->ArcCenter(), center2 = arc2->ArcCenter();
			double l = distance(center1, center2);
			if (abs(l) < 4 * JudgeACC) return false;
			dir = Point2d((center1.x + center2.x) / 2, (center1.y + center2.y) / 2);
			return true;
		}

		void removering(Ring* r) {
			auto it = std::find_if(rings.begin(), rings.end(), [r](const Ring* element) {return element == r;});
			if (it != rings.end()) {
				Ringfordiscard.push_back(*it); rings.erase(it);
			}
		}
		void changebase(std::vector<Position> p, Vertex* oldbase, Vertex* newbase) {
			if (oldbase == newbase) return;
			for (auto& pp : p) {
				if (pp.base == oldbase) pp.base = newbase;
			}
		}
		IntervalsforL getL() {
			auto lastd2 = -INFINITY; std::deque <Point2d> inters; std::vector<int> Beremove;
			std::vector<Vertex*> interself;
			std::sort(lefts.begin(), lefts.end(), [](const Item& A, const Item& B) {return A.dist < B.dist;});
			std::sort(rights.begin(), rights.end(), [](const Item& A, const Item& B) {return A.dist < B.dist;});
			for (int i = 0;i < fmin(lefts.size(), rights.size());++i) {
				double d1 = lefts[i].dist, d2 = rights[i].dist;
				if (lefts[i].dist < lastd2) {//self-intersect
					if (lefts[i].dist + radius < JudgeACC && rights[i].dist - radius > JudgeACC) { //current zero
						Beremove.push_back(i - 1);
						inters.pop_back();//discard last interval;
					}
					else {
						Beremove.push_back(i);
						continue;//discard this
					}
				}
				inters.emplace_back(lefts[i].dist, rights[i].dist); lastd2 = d2;
			}
			std::vector<double> result;
			for (auto& i : inters) {
				double x = i.x + radius - FloatACC;
				double y = i.y - radius + FloatACC;
				if (x < y) {
					result.push_back(x);result.push_back(y);
				}
			}
			if (Beremove.size()) {
				std::reverse(Beremove.begin(), Beremove.end());
				for (int a : Beremove) {
					interself.push_back(lefts[a].base);
					interself.push_back(rights[a].base);
				}
				for (auto it = Beremove.begin(); it != Beremove.end(); ++it) {
					lefts.erase(lefts.begin() + *it);
					rights.erase(rights.begin() + *it);
				}
			}
			for (auto& v : interself) {
				remove_element(VaildSides, v);
				if (v->concave) remove_element(VaildVertices, v);
				if (v->Below->concave) remove_element(VaildVertices, v->Below);
			}
			return IntervalsforL(result);
		}
		struct Item {
			Item(double dist, Vertex* base, Point2d pt) :dist(dist), base(base), UsePos(true), pt(pt) {}
			Item(double dist, Vertex* base) :dist(dist), base(base), UsePos(false), pt(0, 0) {}
			double dist;Vertex* base;bool UsePos; Point2d pt;
		};
		void addInterval(Vertex* v) {
			if (v->m == 0) {
				addInterval(v->m > v->Below->m, v->n, v, Point2d(v->x, v->y));
			}
			else if (v->Below->m == 0) {
				addInterval(v->m > v->Below->m, v->n, v->Below, Point2d(v->Below->x, v->Below->y));
			}
			if (v->m * v->Below->m < 0) {
				if (v->arc) {
					double R = (v->arc >= CurrentArcNo) ? radius : maximumRadius;
					double dt = R - v->length / 2 < 0 ? 0 : sqrt(R * R - v->length * v->length / 4);
					double cm = (v->m + v->Below->m) / 2 + dt / v->length * (v->Below->n - v->n);
					double cn = (v->n + v->Below->n) / 2 - dt / v->length * (v->Below->m - v->m);
					double t = sqrt(R * R - cm * cm);
					if (v->m > v->Below->m) {
						addInterval(true, cn - t, v, Coor(0, cn - t));
					}
					else {
						addInterval(false, cn + t, v, Coor(0, cn + t));
					}
				}
				else if (fabs(v->n - v->Below->n) < 2 * JudgeACC) {
					addInterval(v->m > v->Below->m, v->n, v, Coor(0, v->n));
				}
			}
		}
		void addInterval(bool left, double dist, Vertex* base, Point2d pt) {
			for (auto& t : left ? lefts : rights) {
				if (fabs(t.dist - dist) < 2 * JudgeACC) return;
			}
			(left ? lefts : rights).emplace_back(dist, base, pt);
		}
		void addInterval(bool left, double dist, Vertex* base) {
			for (auto& t : left ? lefts : rights) {
				if (fabs(t.dist - dist) < 2 * JudgeACC) {
					if (t.UsePos) {
						t.UsePos = false;t.base = base;
					}
					return;
				}
			}
			(left ? lefts : rights).emplace_back(dist, base);
		}
		std::vector<Item> lefts, rights;
		IntervalsforL InitFilter() {
			VaildVertices.clear();VaildSides.clear();lefts.clear(); rights.clear();
			for (Ring* r : rings) {
				auto f = r->first;
				project(f);project(f->Below);
				if (f->m >= 0 == f->Below->m >= 0) {
					searchotherring(r);
				}
				else {
					searchfirstring(r);  //intersection on first side
				}
			}
			return getL();
		}

		std::vector<Position> getPreciseStadiumContacts(IntervalsforL& pp, double& testL, bool& isNormal) {
			auto Max = pp.Max(); testL = (Max.x + Max.y) / 2; 	isNormal = false;
			auto p = getStadiumContacts(Max.x, Max.y); isNormal = Max.y - Max.x < 2 * JudgeACC;

			if (p.size() == 1) return p;
			for (int i = 0; i < p.size();i++) {
				Position f = p[i]; int tag1 = distance(*f.base, f) < JudgeACC ? 0 : f.base->arc + 1;
				for (int j = i + 1;j < p.size();j++) {
					Position t = p[j]; int r = 0;
					if (f.base == t.base) continue;
					int tag2 = distance(*t.base, t) < JudgeACC ? 0 : t.base->arc + 1;
					if (tag1 == 0) {
						if (tag2 == 0) {
							r = PointPoint(testL, radius, f.base, t.base);
						}
						else if (tag2 == 1) {
							r = PointLine(testL, radius, f.base, t.base);
						}
						else {
							r = PointCircle(testL, radius, f.base, t.base);
						}
					}
					else if (tag1 == 1) {
						if (tag2 == 0) {
							r = PointLine(testL, radius, t.base, f.base);
						}
						else if (tag2 == 1) {
							r = LineLine(testL, radius, f.base, t.base);
						}
						else {
							r = LineCircle(testL, radius, f.base, t.base);
						}
					}
					else {
						if (tag2 == 0) {
							r = PointCircle(testL, radius, t.base, f.base);
						}
						else if (tag2 == 1) {
							r = LineCircle(testL, radius, t.base, f.base);
						}
						else {
							r = CircleCircle(testL, radius, t.base, f.base);
						}
					}
					if (r == 1) {
						isNormal = true;
						p = getStadiumContacts(testL, testL);
					}
				}
			}
			return p;
		}
		std::vector<Position> getStadiumContacts(double n1, double n2) {
			std::vector<Position> Contacts;
			double radiusJudge = radius + 4 * JudgeACC;
			for (const auto& v : VaildVertices) {
				if ((v->n - n1 - JudgeACC) * (n2 - v->n + JudgeACC) > 0) {
					if (v->m > 0 && v->m < radiusJudge) {   //up
						Position p(v);p.note = (n2 - v->n) / (n2 - n1);
						Contacts.push_back(p);
					}
					else if (v->m < 0 && -v->m < radiusJudge) { //down
						Position p(v);p.note = 2 + (v->n - n1) / (n2 - n1);
						Contacts.push_back(p);
					}
				}
				else if (v->n < n1 + JudgeACC && sqrt((v->n - n1) * (v->n - n1) + v->m * v->m) < radiusJudge) {
					Position p(v);p.note = 1.5 - v->m / radius / 2;
					Contacts.push_back(p);
				}
				else if (v->n > n2 - JudgeACC && sqrt((v->n - n2) * (v->n - n2) + v->m * v->m) < radiusJudge) {
					Position p(v);p.note = 3.5 + v->m / radius / 2;
					if (p.note > 4 - JudgeACC) p.note = 0;
					Contacts.push_back(p);
				}

			}
			for (const auto& v : VaildSides) {
				if (v->arc) {
					double cm = v->d, cn = v->sinA;
					if ((v->d1 - n1 - JudgeACC) * (n1 - v->d2 + JudgeACC) > 0 && v->m >= v->Below->m) {
						double dist = sqrt(cm * cm + (cn - n1) * (cn - n1));
						if (maximumRadius - dist < radiusJudge) {
							auto pt = Coor(-radius * cm / dist, n1 + radius * (n1 - cn) / dist);
							if (atLeft(pt, *v->Below, *v, JudgeACC)) {
								Position p(v, pt.x, pt.y); p.note = 1.5 + cm / dist / 2;
								Contacts.push_back(p);
							}
						}
					}
					if ((v->d1 - n2 - JudgeACC) * (n2 - v->d2 + JudgeACC) > 0 && v->m < v->Below->m) {
						double dist = sqrt(cm * cm + (cn - n2) * (cn - n2));
						if (maximumRadius - dist < radiusJudge) {
							auto pt = Coor(-radius * cm / dist, n2 + radius * (n2 - cn) / dist);
							if (atLeft(pt, *v->Below, *v, JudgeACC)) {
								Position p(v, pt.x, pt.y); p.note = 3.5 - cm / dist / 2;
								if (p.note > 4 - JudgeACC) p.note = 0;
								Contacts.push_back(p);
							}
						}
					}
				}
				else if (n2 - n1 > 2 * JudgeACC && fabs((n2 - n1) * v->sinA / v->cosA) < JudgeACC) {
					if ((v->d1 - n1 - JudgeACC) * (n1 - v->d2 + JudgeACC) > 0) {
						auto pt = Coor(v->m + (n1 - v->n) * v->sinA / v->cosA, n1);
						Position p(v, pt.x, pt.y); p.note = v->cosA < 0 ? 1 : 2;
						Contacts.push_back(p);
					}
					if ((v->d1 - n2 - JudgeACC) * (n2 - v->d2 + JudgeACC) > 0) {
						auto pt = Coor(v->m + (n2 - v->n) * v->sinA / v->cosA, n2);
						Position p(v, pt.x, pt.y); p.note = v->cosA < 0 ? 0 : 3;
						Contacts.push_back(p);
					}
					double n = (n1 + n2) / 2;
					if ((v->d1 - n - JudgeACC) * (n - v->d2 + JudgeACC) > 0) {
						auto pt = Coor(v->m + (n - v->n) * v->sinA / v->cosA, n);
						Position p(v, pt.x, pt.y); p.note = v->cosA < 0 ? 0.5 : 2.5;
						Contacts.push_back(p);
					}
				}
				else if (v->sinA < 0) { //left
					if ((v->d1 - n1 - JudgeACC) * (n1 + JudgeACC - v->d2) < 0) continue;
					double r = fabs(v->d + n1 * v->sinA);
					if (r < radiusJudge) {
						double l = ((n1 - v->n) * (v->Below->n - v->n) - v->m * (v->Below->m - v->m)) / v->length;
						Position p(v, v->x + l * v->cosL, v->y + l * v->sinL); p.note = 1.5 + v->cosA / 2;
						Contacts.push_back(p);
					}
				}
				else {
					if ((v->d1 - n2 - JudgeACC) * (n2 + JudgeACC - v->d2) < 0) continue;
					double r = fabs(v->d + n2 * v->sinA);
					if (r < radiusJudge) {
						double l = ((n2 - v->n) * (v->Below->n - v->n) - v->m * (v->Below->m - v->m)) / v->length;
						Position p(v, v->x + l * v->cosL, v->y + l * v->sinL); p.note = 3.5 - v->cosA / 2;
						if (p.note > 4 - JudgeACC) p.note = 0;
						Contacts.push_back(p);
					}
				}
			}

			addInitContacts(Contacts, n1, n2);
			return Contacts;
		}
		void addInitContacts(std::vector<Position>& result, double n1, double n2) {
			double acc = fmax(JudgeACC, sqrt(2 * maximumRadius * JudgeACC));
			for (auto& a : lefts) {
				if (a.UsePos && fabs(a.dist + radius - n1) < acc) {
					Position p(a.base, a.pt.x, a.pt.y); p.note = 1.5;
					result.push_back(p);
				}
			}
			for (auto& a : rights) {
				if (a.UsePos && fabs(a.dist - radius - n2) < acc) {
					Position p(a.base, a.pt.x, a.pt.y); p.note = 3.5;
					result.push_back(p);
				}
			}
			if (!result.empty()) {
				std::sort(result.begin(), result.end(), [](const Position& A, const Position& B) { return A.note < B.note;});
				if (result.size() > 1 && distance(result.front(), result.back()) < JudgeACC) result.pop_back();
				for (auto it = result.begin();it < result.end() - 1;) {
					auto nextit = std::next(it);
					if (distance(*it, *nextit) < 2 * JudgeACC) {
						result.erase(nextit);
					}
					else {
						if (distance(*it->base->Below, *it) < JudgeACC) it->base = it->base->Below;
						++it;
					}
				}    //Eliminate duplicates
			}
		}
		Point2d Coor(double m, double n) {
			return Point2d(x0 + n * cosR - m * sinR, y0 + m * cosR + n * sinR);
		}
		std::deque<Ring*>& rings; std::vector<Vertex*>VaildVertices, VaildSides;
		void Constraint(IntervalsforL& pp, double testR) {
			std::vector<Vertex*> Vs, Ls;
			for (auto& v : VaildVertices) {
				if (pp.isVertexUseful(v, testR)) {
					Vs.push_back(v); pp.vertexConstraint(v, testR);
				}
			}
			for (auto& v : VaildSides) {
				if (pp.isSideUseful(v, testR)) {

					Ls.push_back(v); pp.sideConstraint(v, testR);
				}
			}
			VaildVertices = Vs; VaildSides = Ls;
			for (auto& v : lefts) {
				pp.Minus(v.dist - JudgeACC, v.dist + testR);
			}
			for (auto& v : rights) {
				pp.Minus(v.dist - testR, v.dist + JudgeACC);
			}
		}
		double GetRadius(double testL) {
			double R = INFINITY;bool f = false;
			for (const auto& v : VaildVertices) {
				v->testRadiusforVertex(testL, R);
			}
			for (const auto& v : VaildSides) {
				v->testRadiusforSide(testL, R);
			}

			for (auto& v : lefts) {
				if (testL > v.dist && R > testL - v.dist) R = testL - v.dist;
			}
			for (auto& v : rights) {
				if (testL<v.dist && R> v.dist - testL) R = v.dist - testL;
			}
			return R;
		}
		void project(Vertex* v) const {
			v->m = (v->y - y0) * cosR - (v->x - x0) * sinR;
			v->n = (v->x - x0) * cosR + (v->y - y0) * sinR;
			if (fabs(v->n) < JudgeACC) v->n = 0.0;
			if (fabs(v->m) < JudgeACC) v->m = 0.0;
		}
		int isArcVaild(Vertex* f) {
			if (abs(f->n - f->Below->n) < 2 * JudgeACC) return false;
			double m = f->m, n = f->n, m2 = f->Below->m, n2 = f->Below->n;
			double R = (f->arc >= CurrentArcNo) ? radius : maximumRadius;
			double dt = R - f->length / 2 < 0 ? 0 : sqrt(R * R - f->length * f->length / 4);
			double cm = (m + m2) / 2 + dt / f->length * (n2 - n);
			double cn = (n + n2) / 2 - dt / f->length * (m2 - m);
			if (abs(cm) < 2 * JudgeACC) return false;
			f->sinA = cn; f->d = cm;
			if (m * m2 >= 0) {
				if ((m + m2) * cm < 0) {
					f->d1 = n - m * (cn - n) / (cm - m);
					f->d2 = n2 - m2 * (cn - n2) / (cm - m2);
					f->r1 = fabs(R * m / (cm - m));
					f->r2 = fabs(R * m2 / (cm - m2));
					return 1;
				}
				else {
					return 0;
				}
			}
			else if (m > 0) {
				f->d1 = cn - (R - cm < 0 ? 0 : sqrt(R * R - cm * cm));
				f->r1 = 0.0;
				if (cm > 0) {
					f->d2 = n2 - m2 * (cn - n2) / (cm - m2);
					f->r2 = fabs(R * m2 / (cm - m2));
				}
				else {
					f->d2 = n - m * (cn - n) / (cm - m);
					f->r2 = fabs(R * m / (cm - m));
				}
				addInterval(true, f->d1, f);

				return 2;
			}
			else {
				f->d2 = cn + (R - cm < 0 ? 0 : sqrt(R * R - cm * cm));
				f->r2 = 0;
				if (cm < 0) {
					f->d1 = n2 - m2 * (cn - n2) / (cm - m2);
					f->r1 = fabs(R * m2 / (cm - m2));
				}
				else {
					f->d1 = n - m * (cn - n) / (cm - m);
					f->r1 = fabs(R * m / (cm - m));
				}
				addInterval(false, f->d2, f);
				return 2;
			}
		}
		bool isLineVaild(Vertex* f) {
			if (abs(f->n - f->Below->n) < 2 * JudgeACC) return false;
			if (f->m * f->Below->m >= 0) {
				if ((f->m + f->Below->m) * (f->n - f->Below->n) > 0) {
					f->projectforLine(sinR, cosR);
					return true;
				}
				else {
					return false;
				}
			}
			else {
				f->projectforLine(sinR, cosR);
				double tp = (f->n - f->Below->n) / (f->Below->m - f->m) * f->m + f->n;
				if (f->m > 0 == (f->n - f->Below->n) > 0) {
					f->r2 = 0;f->d2 = tp;
				}
				else {
					f->r1 = 0;f->d1 = tp;
				}
				addInterval(f->m > 0, tp, f);
				return true;
			}
		}

		void searchotherring(Ring* r) {
			auto f = r->first; project(f->Above); bool last = true;
			do {
				project(f->Below); 	addInterval(f);
				if (f->arc ? isArcVaild(f) : isLineVaild(f)) {
					VaildSides.push_back(f);last = true;
					if (f->concave) VaildVertices.push_back(f);
				}
				else {
					if (f->concave && last) VaildVertices.push_back(f);
					last = false;
				}
				f = f->Below;
			} while (f != r->first);
		}
		void searchfirstring(Ring* r) {
			auto f = r->first->Below, e = f;// project(f->Above);project(f);
			double isup = f->m >= 0 ? 1 : -1;
			auto fs = searchtoNext(f, e, isup);
			do {
				isup = -isup; auto temp = f;

				auto ts = searchtoNext(f, e, isup);
				if (!ts.empty() && ts.front() == temp) {
					if (!fs.empty() && fs.back()->Below == temp) {//join
						fs.pop_back();ts.pop_front();
					}
					else {
						ts.front() = temp->Above;
					}
				}
				else if (!fs.empty() && fs.back()->Below == temp) {//join
					fs.back() = temp;
				}
				else {
					fs.push_back(temp->Above);
					fs.push_back(temp);
				}
				fs.insert(fs.end(), ts.begin(), ts.end());
			} while (f != e);
			if (fs.front() == f) {
				fs.front() = f->Above;
			}
			else if (fs.back()->Below == f) {
				fs.back() = f;
			}
			else {
				fs.push_back(f->Above);
				fs.push_back(f);
			}
			for (auto it = fs.begin();it != fs.end();) {
				auto from = *it++, to = *it++;
				do {
					if (from->arc) {
						if (isArcVaild(from)) VaildSides.push_back(from);
					}
					else {
						if (isLineVaild(from)) VaildSides.push_back(from);
					}

					if (from->concave) VaildVertices.push_back(from);
					from = from->Below;
				} while (from != to);
				if (from->concave) VaildVertices.push_back(from);
			}
			if (VaildVertices.size() > 1 && VaildVertices.back() == VaildVertices.front()) VaildVertices.pop_back();
		}
		bool isConcave(Vertex* f) {
			return (f->arc + f->concave) && atLeft(*f->Below, *f, *f->Above, JudgeACC);
		}
		class checkvaild {
		public:
			std::deque<Vertex*> List; Vertex* begin = nullptr, * end = nullptr;bool currentup = true, atleft = true;
			checkvaild(double isup) :isup(isup) {}
			void add(Vertex* newend) {
				if (newend != begin) { List.push_back(begin);List.push_back(newend);lefts.push(atleft);ups.push(currentup); }
			}
			void stopback(double tp) {
				while ((end->Above->n - tp) * isup < JudgeACC && end != begin) {
					end = end->Above;
				}
				add(end);
			}
			bool backward(double ma, double na, double m, double n) {
				while (!List.empty() && atleft == lefts.top()) {
					if (currentup && (n - List.back()->n) * isup > -JudgeACC) {
						auto bk = List.back();
						auto mm = (ma - m) / (na - n) * (bk->n - n) + m;
						if ((mm - bk->m) * isup < 0) {
							end = bk;List.pop_back();
							begin = List.back();List.pop_back();
							currentup = ups.top(); ups.pop();
							atleft = lefts.top();lefts.pop();
						}
						else {
							return true;
						}
					}
					else if (!currentup && (n - begin->n) * isup > -JudgeACC) {
						end = List.back();List.pop_back();
						begin = List.back();List.pop_back();
						currentup = ups.top(); ups.pop();
						atleft = lefts.top();lefts.pop();
					}
					else {
						return false;
					}
				}
				return false;
			}
		private:
			std::stack<bool> ups;std::stack<bool> lefts;double isup;
		};
		std::deque<Vertex*> searchtoNext(Vertex*& f, Vertex* e, double& isup) {
			addInterval(f->Above);
			double temp = 0, left = 0, right = 0; auto cv = checkvaild(isup);
			int forward = 0;	//status£º 1 advance£¬2 recoil£¬ 3 up£¬4 uptoright £¬5 crinkle
			auto intN = f->Above->arc ? f->n : (f->Above->n - f->n) / (f->m - f->Above->m) * f->Above->m + f->Above->n;
			left = intN;right = intN;
			if (f->Above->arc) {
				f = f->Below;   project(f);
				if (isup > 0 != f->m >= 0)return cv.List;//The first side crosses the line;
				if ((f->n - f->Above->n) * isup < JudgeACC) {
					forward = 1;cv.begin = f->Above;
				}
				else {
					forward = 4;
				}

			}
			else if ((f->n - intN) * isup < -JudgeACC) {
				forward = 1;cv.begin = f;
			}
			else {
				forward = 4;
			}
			do {
				f = f->Below;  project(f);
				if (isup > 0 != f->m >= 0) {
					if (forward == 2) {
						auto tp = f->Above->arc ? f->Above->n : (f->Above->n - f->n) / (f->m - f->Above->m) * f->Above->m + f->Above->n;
						cv.backward(0, tp - isup, 0, tp);
						tp = isup * (f->Above->n - tp) > 0 ? f->Above->n : tp; //returned position; Arc may discard for self-intersect,so its intersection cannot be used.
						cv.stopback(tp);
					}
					else if (forward == 1) {
						cv.add(f->Above);
						if (cv.List.empty()) {
							auto v = f->Above;
							if (v->arc ? isArcVaild(v) : isLineVaild(v)) VaildSides.push_back(v);
						}
					}

					return cv.List;
				}
				switch (forward)
				{
				case 1: {
					if ((f->n - f->Above->n) * isup > JudgeACC) {
						if (isConcave(f->Above)) {
							cv.add(f->Above);
							if ((f->Above->n - left) * isup < JudgeACC) {
								forward = 4; left = f->Above->n;
							}
							else {//1-3
								forward = 3; temp = f->Above->n;
							}
						}
						else {
							if ((f->Above->n - left) * isup < -JudgeACC) left = f->Above->n;
							cv.end = f->Above; forward = 2; goto f2;
						}
					}
					break;
				}
				case 2: {
					if ((f->n - f->Above->n) * isup < -JudgeACC) {
						if (isConcave(f->Above)) { //2-1
							cv.stopback(f->Above->n); cv.currentup = false;
							cv.begin = f->Above; forward = 1;
						}
						else {//2-5
							forward = 5;temp = f->Above->n;
						}

					}
					else {
					f2:
						if ((f->n - right) * isup > JudgeACC) { forward = 4;break; }
						if (cv.backward(f->Above->m, f->Above->n, f->m, f->n)) { forward = 3; temp = cv.List.back()->n; break; }
					}
					break;
				}
				case 3: {
					if ((f->n - temp) * isup < -JudgeACC) {
						cv.begin = f->Above; cv.currentup = true; forward = 1;
					}
					break;
				}
				case 4: {
					if ((f->n - f->Above->n) * isup < -JudgeACC) {
						if ((f->Above->n - right) * isup > JudgeACC) {
							right = f->Above->n;
							if (isConcave(f->Above)) { //right 4-1
								cv.begin = f->Above;cv.atleft = false;cv.currentup = true;
								forward = 1;break;
							}
						}
						if ((f->n - left) * isup < -JudgeACC) {//left 4-1
							cv.begin = f->Above;cv.atleft = true;cv.currentup = true;
							forward = 1;break;
						}
					}
					break;
				}
				case 5: {
					if ((f->n - temp) * isup > JudgeACC) {
						forward = 2;goto f2;
					}
					break;
				}
				}
			} while (f != e);

			return cv.List;
		}
	};
	class Polygon {
	public:
		Polygon(std::vector<std::vector<Vertex*>> pts) {
			CurrentRingNo = 0;
			for (auto& t : pts) {
				auto r = new Ring(t);
				area += r->area.x + r->area.y;
				rings.push_back(r);
			}
			CurrentArcNo = 0;
		}
		Polygon(const std::deque<Ring*>& rings) : rings(rings) {
			for (auto& r : rings) {
				area += r->area.x + r->area.y;
			}
		}
		Circle maximumCircle() {
			maximumRadius = 0; maximumArea = 0; maximumCenterX = 0; maximumCenterY = 0;
			searchMaximumCircle();
			return Circle(maximumCenterX, maximumCenterY, maximumRadius);
		}
	public:
		std::deque<Ring*> rings; double area = 0;
		void searchMaximumCircle() {
			CurrentArcNo += 1;
			auto line = Datumline(rings);
			auto Piece = line.GetPiece();
			Piece.insert(Piece.begin(), rings.begin(), rings.end());
			std::vector<Ring*>  Outs, Holes;
			for (Ring* r : Piece) {
				if (r->area.x + r->area.y > maximumArea) Outs.push_back(r);
				else if (r->area.x < -JudgeACC) Holes.push_back(r);
			}
			std::sort(Outs.begin(), Outs.end(), [](const Ring* A, const Ring* B) {return (A->area.x + A->area.y) * (A->allowRadius == maximumRadius ? 10 : 1) > (B->area.x + B->area.y) * (B->allowRadius == maximumRadius ? 10 : 1);});
			std::sort(Holes.begin(), Holes.end(), [](const Ring* A, const Ring* B) {return A->area.x + A->area.y < B->area.x + B->area.y;});


			for (Ring* r : Outs) {
				for (Ring* rs : trim(r, Holes)) {
					if (rs->area.x + rs->area.y > maximumArea) {
						auto pg = Polygon(rs, Holes);
						if (pg.area > maximumArea)pg.searchMaximumCircle();
					}
				}
			}
		}
		std::vector<Ring*> trim(Ring* r, std::vector<Ring*>& Holes) {
			std::vector<Ring*> result;

			if (!r->notrim && r->allowRadius < maximumRadius - JudgeACC) {
				int i = r->getClipper();
				if (i == 0) return result;
				if (i == 2) return r->trim(Holes);
			}
			result.push_back(r); return result;
		}
		Polygon(Ring* main, std::vector<Ring*>& Holes) {
			main->isOut = true;
			rings.push_back(main); area = main->area.x + main->area.y;
			for (auto hh = Holes.begin();hh < Holes.end();) {
				auto h = *hh;
				if (-h->area.x - h->area.y < area && main->Contains(h->first)) {
					area -= h->area.x + h->area.y;
					rings.push_back(h); hh = Holes.erase(hh);
				}
				else {
					hh++;
				}
			}
		}
	};
}