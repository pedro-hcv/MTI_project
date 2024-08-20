# include "tool.h"


static glm::dvec2 solve_quad(double a, double b, double c){
    double x1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
    double x2 = (-b - sqrt(b*b - 4*a*c))/(2*a);
    
    return glm::dvec2(x1, x2);
}

static double myatan2(double y, double x) {
	double t = atan2(y,x);
	if (t > 0.) {
		return t;
	} else {
		return t + 2*M_PI;
	}
}

glm::dvec2 tool::line::closest_point(glm::dvec2 xq) const {

	double t_location = ((xq.x - initial_vector.x)*parallel_vector.x + (xq.y - initial_vector.y)*parallel_vector.y)/(parallel_vector.x*parallel_vector.x + parallel_vector.y*parallel_vector.y);
	double px = parallel_vector.x*t_location + initial_vector.x; 
	double py = parallel_vector.y*t_location + initial_vector.y; 
	return glm::dvec2(px, py);
}


// TODO intersect lines:
glm::dvec2 tool::line::intersect(tool::line l) const {

	// parallel vertical lines
	if (abs(parallel_vector.x*l.parallel_vector.y - parallel_vector.y*l.parallel_vector.x) < 1e-12){
		return glm::dvec2(DBL_MAX, DBL_MAX);
	}
	else if (parallel_vector.x == 0){
		double s = (initial_vector.x - l.initial_vector.x)/parallel_vector.x;
		double x = (s*l.parallel_vector.x) + l.initial_vector.x;
		double y = (s*l.parallel_vector.y) + l.initial_vector.y;
		return glm::dvec2(x, y);
	}
	else{
		double s = (l.initial_vector.y - initial_vector.y - (parallel_vector.y/parallel_vector.x)*(l.initial_vector.x - initial_vector.x))/((parallel_vector.y/parallel_vector.x)*l.parallel_vector.x - l.parallel_vector.y);
		double x = (s*l.parallel_vector.x) + l.initial_vector.x;
		double y = (s*l.parallel_vector.y) + l.initial_vector.y;
		return glm::dvec2(x, y);
	}

}

double tool::line::angle_between_lines(tool::line l) const{
	return acos(glm::dot(parallel_vector,l.parallel_vector)/(glm::length(parallel_vector)*glm::length(l.parallel_vector)));;
}


double tool::line::get_length() const{
	return length;
}

tool::line::line() {}

tool::line::line(glm::dvec2 p1, glm::dvec2 p2) {
	parallel_vector.x = p2.x - p1.x;
	parallel_vector.y = p2.y - p1.y;

	initial_vector = p1;

}

tool::line::line(glm::dvec2 p1, glm::dvec2 p2, bool segment):segment(segment) {
	parallel_vector.x = p2.x - p1.x;
	parallel_vector.y = p2.y - p1.y;

	initial_vector = p1;

	if (segment){
		length = glm::length(p2-p1);
		limits = glm::dvec2(0,1);
	}

}

bool tool::bbox::in(glm::dvec2 qp) {
	bool in_x = qp.x >= bbmin_x && qp.x <= bbmax_x;
	bool in_y = qp.y >= bbmin_y && qp.y <= bbmax_y;
	return in_x && in_y;
}

bool tool::bbox::valid() const {
	bool invalid_x = bbmax_x - bbmin_x  < 1e-12;
	bool invalid_y = bbmax_y - bbmin_y  < 1e-12;

	return !(invalid_x || invalid_y);
}

tool::bbox::bbox() {}

tool::bbox::bbox(glm::dvec2 p1, glm::dvec2 p2) {
	bbmin_x = fmin(p1.x, p2.x);
	bbmax_x = fmax(p1.x, p2.x);
	bbmin_y = fmin(p1.y, p2.y);
	bbmax_y = fmax(p1.y, p2.y);
}

tool::bbox::bbox(double bbmin_x, double bbmax_x, double bbmin_y, double bbmax_y) :
				bbmin_x(bbmin_x), bbmax_x(bbmax_x), bbmin_y(bbmin_y), bbmax_y(bbmax_y) {}


tool::circle_segment::circle_segment(double r, double t1, double t2, glm::dvec2 p) : r(r), t1(t1), t2(t2), p(p) {}

tool::circle_segment::circle_segment(glm::dvec2 p1, glm::dvec2 p2, glm::dvec2 p3) {
	double x1 = p1.x; double y1 = p1.y;
	double x2 = p2.x; double y2 = p2.y;
	double x3 = p3.x; double y3 = p3.y;

	glm::dmat3x3 d1(x1*x1 + y1*y1, x2*x2 + y2*y2, x3*x3 + y3*y3, y1, y2, y3, 1. ,1., 1.);
	glm::dmat3x3 d2(x1, x2, x3, x1*x1 + y1*y1, x2*x2 + y2*y2, x3*x3 + y3*y3, 1. ,1., 1.);
	glm::dmat3x3 frac(x1, x2, x3, y1, y2, y3, 1., 1., 1.);

	p.x = glm::determinant(d1)/(2.*glm::determinant(frac));
	p.y = glm::determinant(d2)/(2.*glm::determinant(frac));

	r = glm::length(p-p1);
}

/* unsigned int tool::circle_segment::intersect(glm::dvec2 p1, glm::dvec2 p2, glm::dvec2 &out_i1, glm::dvec2 &out_i2) {
	//c.f. http://stackoverflow.com/questions/1073336/circle-line-segment-collision-detection-algorithm
	// answer by @multitaskPro w/ comments by @Duq

	glm::dvec2 c = p;
	glm::dvec2 i1,i2;

	if (fabs(p1.x - p2.x) > 1e-12) {

		glm::dvec2 p3(p1.x - c.x, p1.y - c.y);
		glm::dvec2 p4(p2.x - c.x, p2.y - c.y);

		double m = (p4.y - p3.y) / (p4.x - p3.x);
		double b = p3.y - m * p3.x;

		double under_radical = r*r * (m*m + 1) - b*b;

		if (under_radical < 0.) return 0;

		double tt1 = (-2*m*b+2*sqrt(under_radical))/(2*m*m + 2);
		double tt2 = (-2*m*b-2*sqrt(under_radical))/(2*m*m + 2);

		i1 = glm::dvec2(tt1 + c.x, m*tt1+b+c.y);
		i2 = glm::dvec2(tt2 + c.x, m*tt2+b+c.y);
	} else {	//vertical case
		double vert_x = p1.x;
		if (vert_x < c.x - r) return 0;
		if (vert_x > c.x + r) return 0;

		// gogo gadget pythagoras
		double dist   = fabs(vert_x-c.x);
		assert(r*r >= dist*dist);
		double height = sqrt(r*r-dist*dist);

		i1 = glm::dvec2(vert_x,  height + c.y);
		i2 = glm::dvec2(vert_x, -height + c.y);
	}

	double ti1 = myatan2(c.y-i1.y, c.x-i1.x);
	bool valid_ti1 = ti1 > fmin(t1, t2) && ti1 < fmax(t1, t2);

	double ti2 = myatan2(c.y-i2.y, c.x-i2.x);
	bool valid_ti2 = ti2 > fmin(t1, t2) && ti2 < fmax(t1, t2);

	if (!valid_ti1 && !valid_ti2) return 0;

	if (valid_ti1 && !valid_ti2) {
		out_i1 = i1;
		return 1;
	}

	if (!valid_ti1 && valid_ti2) {
		out_i1 = i2;
		return 1;
	}

	if (valid_ti1 && valid_ti2) {
		out_i1 = i1;
		out_i2 = i2;
		return 2;
	}

	assert(0);
	return -1;	//unreachable
} */

tool::circle_segment::circle_segment() {}



double tool::circle_segment::distance(glm::dvec2 qp) const {

	double t = myatan2(qp.y - p.y, qp.x - p.x); 

	bool valid = t > fmin(t1,t2) && t < fmax(t1,t2);

	if (valid) {
		return sqrt((p.x + r*(cos(t)) - qp.x)*(p.x + r*(cos(t)) - qp.x) + (p.y + r*(sin(t)) - qp.y)*(p.y + r*(sin(t)) - qp.y));
	} else {
		return DBL_MAX;
	}
	//glm::dvec2 cross();
}

glm::dvec2 tool::fit_fillet(double r, tool::line lm, tool::line l1) const {

	glm::dvec2 intersec = lm.intersect(l1);
	//double gamma = acos();
	// angle between lines:
	double alpha = lm.angle_between_lines(l1);
	double xm = intersec.x - r;
	double ym = intersec.y + (r/tan(alpha/2));
	double gamma = acos((-lm.parallel_vector.y)/(glm::length(lm.parallel_vector)));
	
	if (gamma > 1e-6){
		double dum_x = xm;
		double dum_y = ym;
		xm = dum_x*cos(gamma) - dum_y*sin(gamma);  
		ym = dum_x*sin(gamma) + dum_y*cos(gamma); 
	}
	return glm::dvec2(xm, ym);
}

void tool::construct_segments(std::vector<glm::dvec2> list_p) {
	unsigned int n = list_p.size();
	for (unsigned int i = 0; i < n; i++) {
		unsigned int cur  = i;
		unsigned int next = (cur+1 > n-1) ? 0 : i+1;
		m_segments.push_back(line(list_p[cur], list_p[next]));
		//m_segments[i].left = m
		//		printf("lft: %f %f\n", m_segments[i].left.x, m_segments[i].left.y);
		//		printf("rgt: %f %f\n", m_segments[i].right.x, m_segments[i].right.y);
		//		printf("nrm: %f %f\n", m_segments[i].n.x, m_segments[i].n.y);
		//		printf("ab:  %f %f\n", m_segments[i].l.a, m_segments[i].l.b);
		//		printf("======================\n");
	}
}


std::vector<glm::dvec2> tool::construct_points_and_fillet(glm::dvec2 tl, glm::dvec2 tr, glm::dvec2 br, glm::dvec2 bl, double r) {

	// find center of fillet => p
	line l1(tr, br);
	line l2(bl, br);
	glm::dvec2 p = fit_fillet(r, l1, l2); // fit_fillet(r, lm, l2) would work too

	// find points on l1, l2 that meet the fillet => trc, blc (c = "continued")
	glm::dvec2 trc = l1.closest_point(p);
	glm::dvec2 blc = l2.closest_point(p);


	// construct circle segment
	double t1 = myatan2(trc.y - p.y, trc.x - p.x);
	double t2 = myatan2(blc.y - p.y, blc.x - p.x);
	m_fillet = new circle_segment(r, t1, t2, p);

	return std::vector<glm::dvec2>({tl, tr, trc, blc, bl});
}


double tool::inside(glm::dvec2 qp) const {

	if (m_chamfer_debug) {
		bool in = glm::length(m_fillet->p - qp) < m_fillet->r;
		if (!in) return -1.;
		return m_fillet->distance(qp);
	}

	bool in = false;

	//determine if point is in polygon
	int i, j;
    for (i = 0, j = m_segments.size()-1; i < m_segments.size(); j = i++) {
      if ((((m_segments[i].initial_vector.y <= qp.y) && (qp.y < m_segments[j].initial_vector.y)) ||
    	  ((m_segments[j].initial_vector.y <= qp.y) && (qp.y < m_segments[i].initial_vector.y))) &&
          (qp.x < (m_segments[j].initial_vector.x - m_segments[i].initial_vector.x) * (qp.y - m_segments[i].initial_vector.y) / (m_segments[j].initial_vector.y - m_segments[i].initial_vector.y) + m_segments[i].initial_vector.x))
        in = !in;
    	}
    

	//if fillet is present, qp can also fall into the fillet, so let's check
	if (m_fillet) {
		in = in || glm::length(m_fillet->p - qp) < m_fillet->r;
	}

	//return a dummy value if qp is not inside
	if (!in) return -1.;

	double r = DBL_MAX;

	for (unsigned int i = 0; i < m_segments.size(); i++) {
		double d;
		if (m_fillet && i==2) {					//if there is a fillet present, its always covering the second segment
			d = m_fillet->distance(qp);
		} else {
			glm::dvec2 p = m_segments[i].closest_point(qp);
			d = glm::length(p-qp);
		}

		r = fmin(r, d);
	}

	if (r == 0.) {
		printf("WARNING: point on tool! visibility test unreliable, nudge tool!\n");
	}

	return r;
}


glm::dvec2 tool::project(glm::dvec2 qp) const {
	assert(inside(qp) >= 0.);

	const double nudge = 1e-8;	//TODO: make proportional to dx?
//	const double nudge = 0.0001;

	double min_dist = DBL_MAX;
	glm::dvec2 proj;

	unsigned int seg_number = 0;
	for (auto it = 0; it != m_segments.size(); it++) {

		//TODO
		if (m_chamfer_debug && (seg_number != 2)) {
			seg_number++;
			continue;
		}

		if (seg_number == 2 && m_fillet) {
			double dist = m_fillet->distance(qp);

			if (dist < min_dist) {
				glm::dvec2 dist_vec = glm::normalize(qp-m_fillet->p)*(m_fillet->r+nudge);
				proj = dist_vec + m_fillet->p;
				proj += nudge*glm::normalize(dist_vec);
				min_dist = dist;
			}
		} else {

			glm::dvec2 proj_candidate = m_segments[it].closest_point(qp);
			glm::dvec2 n = proj_candidate - qp;
			double dist = glm::length(n);

			n = n/dist;

			if (dist < min_dist) {
				min_dist = dist;
				proj = proj_candidate + nudge*n;
			}
		}

		seg_number++;
	}

	if (inside(proj) >= 0) {
		printf("projected point inside tool, wtf!\n");
		printf("qp: %f %f\n", qp.x, qp.y);
		printf("mapped to ->\n");
		printf("qp: %f %f\n", proj.x, proj.y);

		printf("---------\n");

		for (auto it = m_segments.begin(); it != m_segments.end(); ++it) {
			printf("%f %f\n", it->left.x, it->left.y);
		}
		printf("%f %f\n", m_segments[0].left.x, m_segments[0].left.y);
	}


	assert(inside(proj) < 0.);	//make sure projection is outside
	return proj;
}

void tool::get_chamfer_data(glm::dvec2 &p, double &r) const {
	p.x = m_fillet->p.x;
	p.y = m_fillet->p.y;
	r   = m_fillet->r;
}

//////////////////////////////////////////////////////
tool::bbox tool::safe_bb(double safety) const {

	double x_min = DBL_MAX; double x_max = -DBL_MAX;
	double y_min = DBL_MAX; double y_max = -DBL_MAX;

	/* for (auto it = m_segments.begin(); it != m_segments.end(); ++it) {
		x_min = fmin(it->left.x, x_min);
		x_max = fmax(it->left.x, x_max);
		y_min = fmin(it->left.y, y_min);
		y_max = fmax(it->left.y, y_max);
	} */
	for (auto it = 0; it != m_segments.size(); it++) {
		x_min = fmin(m_segments[it].initial_vector.x, x_min);
		x_max = fmax(m_segments[it].initial_vector.x, x_max);
		y_min = fmin(m_segments[it].initial_vector.y, y_min);
		y_max = fmax(m_segments[it].initial_vector.y, y_max);
	}

	if (m_fillet) {
		y_min = low();
		x_max = m_fillet->p.x + m_fillet->r;
	}

	x_min -= safety;
	y_min -= safety;
	x_max += safety;
	y_max += safety;

	tool::bbox ret;
	ret.bbmin_x = x_min;
	ret.bbmax_x = x_max;
	ret.bbmin_y = y_min;
	ret.bbmax_y = y_max;

	return ret;
}

double tool::low() const {
	if (m_fillet) {
		return m_fillet->p.y - m_fillet->r;
	}

	double y_min = DBL_MAX;
	for (auto it = m_segments.begin(); it != m_segments.end(); ++it) {
		y_min = fmin(y_min, it->left.y);
	}
	return y_min;
}

//////////////////////////////////////////////////////

bool tool::contact(glm::dvec2 qp, double &depth, glm::dvec2 &dir) const {
	bool in = false;

	//determine if point is in polygon
	int i, j;
    for (i = 0, j = m_segments.size()-1; i < m_segments.size(); j = i++) {
      if ((((m_segments[i].initial_vector.y <= qp.y) && (qp.y < m_segments[j].initial_vector.y)) ||
    	  ((m_segments[j].initial_vector.y <= qp.y) && (qp.y < m_segments[i].initial_vector.y))) &&
          (qp.x < (m_segments[j].initial_vector.x - m_segments[i].initial_vector.x) * (qp.y - m_segments[i].initial_vector.y) / (m_segments[j].initial_vector.y - m_segments[i].initial_vector.y) + m_segments[i].initial_vector.x))
        in = !in;
    	}
    

	//if fillet is present, qp can also fall into the fillet, so let's check
	if (m_fillet) {
		in = in || glm::length(m_fillet->p - qp) < m_fillet->r;
	}


	if (!in) return false;

	depth = DBL_MAX;

	for (unsigned int i = 0; i < m_segments.size(); i++) {
		double d;
		if (m_fillet && i==2) {					//if there is a fillet present, its always covering the second segment
			d = m_fillet->distance(qp);
			if (d < depth) {
				depth = d;
				dir = glm::normalize(qp - m_fillet->p);
			}
		} else {
			glm::dvec2 p = m_segments[i].closest_point(qp);
			d = glm::length(p-qp);
			if (d < depth) {
				depth = d;
				dir = glm::normalize(qp - p);
			}
		}
	}

	return true;
}

bool tool::contact(glm::dvec2 qp, glm::dvec2 &cp, glm::dvec2 &n) const {
	bool in = false;

	//determine if point is in polygon
	int i, j;
    for (i = 0, j = m_segments.size()-1; i < m_segments.size(); j = i++) {
      if ((((m_segments[i].initial_vector.y <= qp.y) && (qp.y < m_segments[j].initial_vector.y)) ||
    	  ((m_segments[j].initial_vector.y <= qp.y) && (qp.y < m_segments[i].initial_vector.y))) &&
          (qp.x < (m_segments[j].initial_vector.x - m_segments[i].initial_vector.x) * (qp.y - m_segments[i].initial_vector.y) / (m_segments[j].initial_vector.y - m_segments[i].initial_vector.y) + m_segments[i].initial_vector.x))
        in = !in;
    	}
    

	//if fillet is present, qp can also fall into the fillet, so let's check
	if (m_fillet) {
		in = in || glm::length(m_fillet->p - qp) < m_fillet->r;
	}


	if (!in) return false;

	double depth = DBL_MAX;

	for (unsigned int i = 0; i < m_segments.size(); i++) {

		if (m_chamfer_debug && i != 2) continue;

		double d;
		if (m_fillet && i==2) {					//if there is a fillet present, its always covering the third segment
			d = m_fillet->distance(qp);
			if (d < depth) {
				depth = d;

				n  = glm::normalize(qp - m_fillet->p);
				cp = m_fillet->p + n*m_fillet->r;
			}
		} else {
			glm::dvec2 p = m_segments[i].closest_point(qp);
			d = glm::length(p-qp);
			if (d < depth) {
				depth = d;

				n  = glm::normalize(qp - p);
				cp = p;
			}
		}
	}

	return true;
}

void tool::update_tool() {
	simulation_time *time = &simulation_time::getInstance();
	double delta_t = time->get_dt();
	update_tool(delta_t);

}

void tool::update_tool(double dt) {
	for (auto it = 0; it != m_segments.size(); it++) {
		//it->left  += dt*m_velocity;
		//it->right += dt*m_velocity;
		//it->l = line(it->left, it->right);
		m_segments[it].initial_vector +=  dt*m_velocity;
	}

	if (m_fillet) {
		m_fillet->p += dt*m_velocity;
	}

	for (auto it = m_boundary_particles.begin(); it != m_boundary_particles.end(); ++it) {
		it->x += dt*m_velocity.x;
		it->y += dt*m_velocity.y;
	}
}

void tool::set_vel(glm::dvec2 vel) {
	m_velocity = vel;
}

glm::dvec2 tool::get_vel() const{
	return m_velocity;
}

void tool::set_edge_coord(glm::dvec2 coord) {
	m_edge_coord = coord;
}

glm::dvec2 tool::get_edge_coord() const{
	return m_edge_coord;
}

glm::dvec2 tool::center() const {
	double cx = 0.;
	double cy = 0.;

	for (auto it = m_segments.begin(); it != m_segments.end(); ++it) {
		cx += it->initial_vector.x;
		cy += it->initial_vector.y;
	}

	cx /= m_segments.size();
	cy /= m_segments.size();

	return glm::dvec2(cx, cy);
}

void tool::print(FILE *fp) {
	for (unsigned int it = 0; it < m_segments.size(); ++it) {
		fprintf(fp, "%f %f\n", m_segments[it].initial_vector.x, m_segments[it].initial_vector.y);
	}
	fprintf(fp, "%f %f\n", m_segments.back().initial_vector.x, m_segments.back().initial_vector.y);

	fprintf(fp, "%f %f\n", m_fillet->p.x, m_fillet->p.y);
	fprintf(fp, "%f %f\n", m_fillet->r, m_fillet->r);
}

void tool::print(unsigned int step, const char *folder_name) {
	char buf[256];
	sprintf(buf, "./%s/tool_%07d.txt", folder_name, step);

	FILE *fp = fopen(buf, "w+");
	if (!fp) {
		printf("could not open tool file\n!");
		exit(-1);
	}

	fprintf(fp, "%d\n", (int) m_segments.size());
	for (auto it = m_segments.begin(); it != m_segments.end(); ++it) {
		fprintf(fp, "%f %f\n", it->initial_vector.x, it->initial_vector.y);
	}

	if (m_fillet) {
		fprintf(fp, "%f %f %f %f %f\n", m_fillet->p.x, m_fillet->p.y, m_fillet->r, m_fillet->t1, m_fillet->t2);
	}

	fclose(fp);
}

double tool::mu() const {
	return m_mu;
}

void tool::set_chamfer(glm::dvec2 cp, double r, double t1, double t2) {
	m_fillet = new circle_segment(r, t1, t2, cp);
}

void tool::set_chamfer_debug(bool chamfer_debug) {
	m_chamfer_debug = chamfer_debug;
}

tool::tool(glm::dvec2 tl, glm::dvec2 tr, glm::dvec2 br, glm::dvec2 bl, double mu_fric) : m_mu(mu_fric) {
	construct_segments(std::vector<glm::dvec2>({tl, tr, br, bl}));
}

tool::tool(glm::dvec2 tl, glm::dvec2 tr, glm::dvec2 br, glm::dvec2 bl, double r, double mu_fric) : m_mu(mu_fric)  {
	if (r == 0.) {
		construct_segments(std::vector<glm::dvec2>({tl, tr, br, bl}));
		return;
	}

	std::vector<glm::dvec2> points = construct_points_and_fillet(tl,tr,br,bl,r);
	construct_segments(points);
}

tool::tool(glm::dvec2 tl, double length, double height,
		double rake_angle, double clearance_angle,
		double r, double mu_fric) : m_mu(mu_fric)  {

	glm::dvec2 tr(tl.x+length, tl.y);
	glm::dvec2 bl(tl.x, tl.y-height);

	double alpha_rake = rake_angle * M_PI / 180.;
	//double alpha_free = (180-90-clearance_angle) * M_PI / 180.;
	double alpha_free = (clearance_angle) * M_PI / 180.;

	/* glm::dmat2x2 rot_rake(cos(alpha_rake), -sin(alpha_rake), sin(alpha_rake), cos(alpha_rake));
	glm::dmat2x2 rot_free(cos(alpha_free), -sin(alpha_free), sin(alpha_free), cos(alpha_free));

	glm::dvec2 down(0., -1.);

	glm::dvec2 trc = tr + down*rot_rake;
	glm::dvec2 blc = bl + down*rot_free;

	tool::line l1(tr,trc);
	tool::line l2(bl,blc);

	glm::dvec2 br = l1.intersect(l2); */
	glm::dvec2 br;
	br.x = tr.x + tan(alpha_rake)*height;
	br.y = bl.y - tan(alpha_free)*length;

	if (r == 0.) {
		construct_segments(std::vector<glm::dvec2>({tl, tr, br, bl}));
		return;
	}

	std::vector<glm::dvec2> points = construct_points_and_fillet(tl,tr,br,bl,r);

	construct_segments(points);
}

const std::vector<tool::line> &tool::get_segments() const {
	return m_segments;
}

const tool::circle_segment *tool::get_fillet() const {
	return m_fillet;
}


tool::tool(glm::dvec2 tl, double length, double height,
		double rake_angle, double clearance_angle,
		double mu_fric) : m_mu(mu_fric)  {
	
	glm::dvec2 tr(tl.x+length, tl.y);
	glm::dvec2 bl(tl.x, tl.y-height);

	double alpha_rake = rake_angle * M_PI / 180.;
	double alpha_free = (clearance_angle) * M_PI / 180.;

	glm::dvec2 br;
	br.x = tr.x + tan(alpha_rake)*height;
	br.y = bl.y - tan(alpha_free)*length;


	construct_segments(std::vector<glm::dvec2>({tl, tr, br, bl}));
}

tool::tool() {}
