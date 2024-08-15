
#ifndef TOOL_H_

#define TOOL_H_ 

#include <glm/glm.hpp>
#include <vector>
#include <stdio.h>
#include <assert.h>


#include "kernel.h"
#include "simulation_time.h"

class tool {
    private:
        struct line
        {
            glm::dvec2 parallel_vector;
            glm::dvec2 initial_vector;
            
            bool segment = false;

            glm::dvec2 limits;

            double length = 0;

            glm::dvec2 closest_point(glm::dvec2 xq) const;

            glm::dvec2 intersect(line l ) const;
            
            glm::dvec2 left;
            
            glm::dvec2 right;

            double angle_between_lines(line l) const;


            double get_length() const;
            
            line(glm::dvec2 pq, glm::dvec2 p2);
            line(glm::dvec2 p1, glm::dvec2 p2, bool segment);

            ///line(double a, double b, bool vertical);
            line();


        };
        

        struct circle_segment
        {
            double r = 0;
            double t1 = 0;
            double t2 = 0;
            glm::dvec2 p;
            
            double distance(glm::dvec2 qp) const;

            unsigned int intersect(glm::dvec2 p1, glm::dvec2 p2, glm::dvec2 &i1, glm::dvec2 &i2);

            circle_segment(double r, double t1, double t2, glm::dvec2 p);
            circle_segment(glm::dvec2 p1, glm::dvec2 p2, glm::dvec2 p3);
            circle_segment();

        };


    private:
        glm::dvec2 fit_fillet(double r, line lm, line l1) const;
        
        void construct_segments(std::vector<glm::dvec2> list_p);

        std::vector<glm::dvec2> construct_points_and_fillet(glm::dvec2 tl, glm::dvec2 tr, glm::dvec2 br, glm::dvec2 bl, double r);

        glm::dvec2 m_velocity;

        glm::dvec2 m_edge_coord;

        double m_mu =0.;

        circle_segment *m_fillet = 0;
        std::vector<line> m_segments;

        std::vector<glm::dvec2> m_boundary_particles;
        std::vector<glm::dvec2> m_boundary_normals;
        std::vector<glm::dvec2> m_boundary_weights;

        bool m_chamfer_debug = false;

    public:
        struct bbox
        {
            double bbmin_x = 0.;
            double bbmax_x = 0.;            
            double bbmin_y = 0.;
            double bbmax_y = 0.;
            
            bool in(glm::dvec2 qp);

            bool valid() const;

            bbox();
            bbox(glm::dvec2 p1, glm::dvec2 p2);
            bbox(double bbmin_x, double bbmax_x, double bbmin_y, double bbmax_y);
            

        };
        
        const std::vector<line> &get_segments() const;
        const circle_segment *get_fillet() const;

        glm::dvec2 project(glm::dvec2 qp) const;

        tool::bbox safe_bb(double safety = 0.011) const;

        double low() const;

        double mu() const;

        double inside(glm::dvec2 qp) const;

        bool intersect(glm::dvec2 p1, glm::dvec2 p2) const;
        bool intersect(glm::dvec2 p1, glm::dvec2 p2, double &r) const;

        bool contact(glm::dvec2 qp, double &depth, glm::dvec2 &dir) const;

        bool contact(glm::dvec2 qp, glm::dvec2 &cp, glm::dvec2 &n) const;

        void update_tool();
        void update_tool(double dt);
        
        void set_vel(glm::dvec2 vel);

        glm::dvec2 get_vel() const;

        void set_edge_coord(glm::dvec2 coord_tip);

        glm::dvec2 get_edge_coord() const;

    	glm::dvec2 center() const;

    	//chamfer debugging
    	void get_chamfer_data(glm::dvec2 &p, double &r) const;
    	void set_chamfer(glm::dvec2 cp, double r, double t1, double t2);
    	void set_chamfer_debug(bool chamfer_debug);

    	// print to file
    	void print(FILE *fp);
    	void print(unsigned int step, const char *folder = "results");

    	// construct tool given by four points and a fillet radius
    	tool(glm::dvec2 tl, glm::dvec2 tr, glm::dvec2 br, glm::dvec2 bl, double r, double mu_fric);
    
    	// construct tool given by four points (tool is perfectly sharp)
    	tool(glm::dvec2 tl, glm::dvec2 tr, glm::dvec2 br, glm::dvec2 bl, double mu_fric);
    
    	// construct tool given by a reference point, length and height
    	// as well as rake and clearance angle (measured from vertically
    	// downwards and horizontally leftwards, respectively), and
    	// fillet radius r
    	// angles are in degrees
    	tool(glm::dvec2 tl, double length, double height,
    			double rake_angle, double clearance_angle,
    			double r, double mu_fric);
    
    	// construct tool given by a reference point, length and height
    	// as well as rake and clearance angle (measured from vertically
    	// downwards and horizontally leftwards, respectively), the tool
    	// is percectly sharp
    	// angles are in degrees
    	tool(glm::dvec2 tl, double length, double height,
    			double rake_angle, double clearance_angle, double mu_fric);
    
    	tool();
};

#endif