namespace iupr_cedges {
struct EdgeDetector {
    virtual void set_gauss(float sx,float sy) = 0;
    virtual void set_noise(float frac,float low,float high) = 0;
    virtual void set_poly(float minlength,float maxdist) = 0;
    virtual void clear() = 0;
    virtual void load_pnm(char *file) = 0;
    virtual void save_pnm(char *file) = 0;
    virtual int dim(int i) = 0;
    virtual void set_image(unsigned char *p,int w,int h) = 0;
    virtual void set_pixmap(unsigned char *p,int w,int h) = 0;
    virtual void compute() = 0;
    virtual void get_eimage(unsigned char *p,int w,int h) = 0;
    virtual void get_epixmap(unsigned char *image,int w,int h) = 0;
    virtual float gradient_magnitude(int x,int y) = 0;
    virtual float gradient_angle(int x,int y) = 0;
    virtual bool nextchain() = 0;
    virtual int npoints() = 0;
    virtual void point(int index,float &x,float &y) = 0;
    virtual int nsegments() = 0;
    virtual void segment(int i,float &x0,float &y0,float &x1,float &y1,
			 float &angle,float &magnitude,int &n) = 0;
    virtual ~EdgeDetector() {}
};

EdgeDetector *makeEdgeDetector();
}
