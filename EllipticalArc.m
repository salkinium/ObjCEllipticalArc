// Copyright (c) 2003-2004, Luc Maisonobe
// All rights reserved.
//
// Redistribution and use in source and binary forms, with
// or without modification, are permitted provided that
// the following conditions are met:
//
//    Redistributions of source code must retain the
//    above copyright notice, this list of conditions and
//    the following disclaimer.
//    Redistributions in binary form must reproduce the
//    above copyright notice, this list of conditions and
//    the following disclaimer in the documentation
//    and/or other materials provided with the
//    distribution.
//    Neither the names of spaceroots.org, spaceroots.com
//    nor the names of their contributors may be used to
//    endorse or promote products derived from this
//    software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
// CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
// THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
// HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
// IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
// USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Ported to Objective-C for CoreGraphics on 12. August 2012 by Niklas Hauser

#import "EllipticalArc.h"

// coefficients for error estimation
// while using quadratic Bézier curves for approximation
// 0 < b/a < 1/4
double coeffs2Low[2][4][4] = {
	{
		{  3.92478,   -13.5822,     -0.233377,    0.0128206   },
		{ -1.08814,     0.859987,    0.000362265, 0.000229036 },
		{ -0.942512,    0.390456,    0.0080909,   0.00723895  },
		{ -0.736228,    0.20998,     0.0129867,   0.0103456   }
	}, {
		{ -0.395018,    6.82464,     0.0995293,   0.0122198   },
		{ -0.545608,    0.0774863,   0.0267327,   0.0132482   },
		{  0.0534754,  -0.0884167,   0.012595,    0.0343396   },
		{  0.209052,   -0.0599987,  -0.00723897,  0.00789976  }
	}
};

// coefficients for error estimation
// while using quadratic Bézier curves for approximation
// 1/4 <= b/a <= 1
double coeffs2High[2][4][4] = {
    {
		{  0.0863805, -11.5595,     -2.68765,     0.181224    },
		{  0.242856,   -1.81073,     1.56876,     1.68544     },
		{  0.233337,   -0.455621,    0.222856,    0.403469    },
		{  0.0612978,  -0.104879,    0.0446799,   0.00867312  }
    }, {
		{  0.028973,    6.68407,     0.171472,    0.0211706   },
		{  0.0307674,  -0.0517815,   0.0216803,  -0.0749348   },
		{ -0.0471179,   0.1288,     -0.0781702,   2.0         },
		{ -0.0309683,   0.0531557,  -0.0227191,   0.0434511   }
    }
};

// safety factor to convert the "best" error approximation
// into a "max bound" error
double safety2[4] = {
    0.02, 2.83, 0.125, 0.01
};

// coefficients for error estimation
// while using cubic Bézier curves for approximation
// 0 < b/a < 1/4
double coeffs3Low[2][4][4] = {
	{
		{  3.85268,   -21.229,      -0.330434,    0.0127842  },
		{ -1.61486,     0.706564,    0.225945,    0.263682   },
		{ -0.910164,    0.388383,    0.00551445,  0.00671814 },
		{ -0.630184,    0.192402,    0.0098871,   0.0102527  }
    }, {
		{ -0.162211,    9.94329,     0.13723,     0.0124084  },
		{ -0.253135,    0.00187735,  0.0230286,   0.01264    },
		{ -0.0695069,  -0.0437594,   0.0120636,   0.0163087  },
		{ -0.0328856,  -0.00926032, -0.00173573,  0.00527385 }
    }
};

// coefficients for error estimation
// while using cubic Bézier curves for approximation
// 1/4 <= b/a <= 1
double coeffs3High[2][4][4] = {
	{
		{  0.0899116, -19.2349,     -4.11711,     0.183362   },
		{  0.138148,   -1.45804,     1.32044,     1.38474    },
		{  0.230903,   -0.450262,    0.219963,    0.414038   },
		{  0.0590565,  -0.101062,    0.0430592,   0.0204699  }
    }, {
		{  0.0164649,   9.89394,     0.0919496,   0.00760802 },
		{  0.0191603,  -0.0322058,   0.0134667,  -0.0825018  },
		{  0.0156192,  -0.017535,    0.00326508, -0.228157   },
		{ -0.0236752,   0.0405821,  -0.0173086,   0.176187   }
    }
};

// safety factor to convert the "best" error approximation
// into a "max bound" error
double safety3[4] = {
    0.001, 4.98, 0.207, 0.0067
};

@interface EllipticalArc()
/** Compute the locations of the focii. */
-(void)computeFocii;

/** Compute the locations of the endpoints. */
-(void)computeEndPoints;

/** Compute the bounding box. */
-(void)computeBounds;

-(void)computeDerivedFlatnessParameters;

/** Compute the value of a rational function.
 * This method handles rational functions where the numerator is
 * quadratic and the denominator is linear
 * @param x abscissa for which the value should be computed
 * @param c coefficients array of the rational function
 */
+(double)rationalFunctionWithAbscissa:(double)x coefficients:(double[])c;

/** Estimate the approximation error for a sub-arc of the instance.
 * @param degree degree of the Bézier curve to use (1, 2 or 3)
 * @param etaA start angle of the sub-arc
 * @param etaB end angle of the sub-arc
 * @return upper bound of the approximation error between the Bézier
 * curve and the real ellipse
 */
-(double)estimateErrorWithDegree:(int)degree startAngle:(double)etaA endAngle:(double)etaB;

/** Tests if a line segment intersects the arc.
 * @param start the first point of the line segment
 * @param end	the second point of the line segment
 * @return true if the two line segments intersect
 */
-(BOOL)intersectLineWithStartPoint:(CGPoint)start endPoint:(CGPoint)end;

/** Tests if two line segments intersect.
 * @param start1	the first point of the first line segment
 * @param end1		the second point of the first line segment
 * @param start2	the first point of the second line segment
 * @param end2		the second point of the second line segment
 * @return true if the two line segments intersect
 */
+(BOOL)lineIntersectLineWithStartPoint1:(CGPoint)start1 endPoint1:(CGPoint)end1
							startPoint2:(CGPoint)start2 endPoint2:(CGPoint)end2;

/** Tests if a line segment intersects the outline.
 * @param start	the first point of the line segment
 * @param end	the second point of the line segment
 * @return true if the two line segments intersect
 */
-(BOOL)outlineIntersectLineWithStartPoint:(CGPoint)start endPoint:(CGPoint)end;

/** Build an approximation of the instance outline.
 * @param degree	degree of the Bézier curve to use
 * @param threshold acceptable error
 * @param m			affine transformation to apply, set to NULL if not desired
 * @return a CGMutablePathRef to be used with CoreGraphics
 */
-(CGMutablePathRef)pathWithDegree:(int)degree threshold:(double)threshold tranformation:(CGAffineTransform *)m;

/** Set the maximal degree allowed for Bézier curve approximation.
 * @param maxDegree maximal allowed degree (must be between 1 and 3)
 */
-(void)setMaxDegree:(int)maxDegree;

/** Set the default flatness for Bézier curve approximation.
 * @param defaultFlatness default flatness (must be greater than 1.0e-10)
 */
-(void)setDefaultFlatness:(double)defaultFlatness;
@end


@implementation EllipticalArc

#pragma mark - constructors
-(id)init
{
	if (self = [super init])
	{
		_cx         = 0;
		_cy         = 0;
		_a          = 1;
		_b          = 1;
		_theta      = 0;
		_eta1       = 0;
		_eta2       = 2 * M_PI;
		_cosTheta   = 1;
		_sinTheta   = 0;
		_isPieSlice = false;
		_maxDegree  = 3;
		_defaultFlatness = 0.5; // half a pixel
		
		[self computeFocii];
		[self computeEndPoints];
		[self computeBounds];
		[self computeDerivedFlatnessParameters];
	}
	return self;
}

-(id)initWithCentre:(CGPoint)centre
	  semiMajorAxis:(double)a
	  semiMinorAxis:(double)b
		orientation:(double)theta
		 startAngle:(double)lambda1
		   endAngle:(double)lambda2
		   pieSlice:(BOOL)pieSlice
{
	if (self = [super init])
	{
		_cx         = centre.x;
		_cy         = centre.y;
		_a          = a;
		_b          = b;
		_theta      = theta;
		_isPieSlice = pieSlice;
		
		_eta1       = atan2(sin(lambda1) / b,
							cos(lambda1) / a);
		_eta2       = atan2(sin(lambda2) / b,
							cos(lambda2) / a);
		_cosTheta   = cos(theta);
		_sinTheta   = sin(theta);
		_maxDegree  = 3;
		_defaultFlatness = 0.5; // half a pixel
		
		// make sure we have eta1 <= eta2 <= eta1 + 2 M_PI
		_eta2 -= (2*M_PI) * floor((_eta2 - _eta1) / (2*M_PI));
		
		// the preceding correction fails if we have exactly et2 - eta1 = 2 M_PI
		// it reduces the interval to zero length
		if ((lambda2 - lambda1 > M_PI) && (_eta2 - _eta1 < M_PI))
		{
			_eta2 += (2*M_PI);
		}
		
		[self computeFocii];
		[self computeEndPoints];
		[self computeBounds];
		[self computeDerivedFlatnessParameters];
	}
	return self;
}

-(id)initWithCentre:(CGPoint)centre
	  semiMajorAxis:(double)a
	  semiMinorAxis:(double)b
		orientation:(double)theta
{
	if (self = [super init])
	{
		_cx         = centre.x;
		_cy         = centre.y;
		_a          = a;
		_b          = b;
		_theta      = theta;
		_isPieSlice = NO;
		
		_eta1      = 0;
		_eta2      = 2 * M_PI;
		_cosTheta  = cos(theta);
		_sinTheta  = sin(theta);
		_maxDegree = 3;
		_defaultFlatness = 0.5; // half a pixel
		
		[self computeFocii];
		[self computeEndPoints];
		[self computeBounds];
		[self computeDerivedFlatnessParameters];
	}
	return self;
}


#pragma mark - setters

-(void)setMaxDegree:(int)maxDegree
{
	if ((maxDegree < 1) || (maxDegree > 3))
		return;
	_maxDegree = maxDegree;
}

-(void)setDefaultFlatness:(double)defaultFlatness
{
	if (defaultFlatness < 1.0e-10)
		return;
	_defaultFlatness = defaultFlatness;
}


#pragma mark - computations

-(void)computeFocii
{
	double d  = sqrt(_a * _a - _b * _b);
    double dx = d * _cosTheta;
    double dy = d * _sinTheta;
	
    _xF1 = _cx - dx;
    _yF1 = _cy - dy;
    _xF2 = _cx + dx;
    _yF2 = _cy + dy;
}

-(void)computeEndPoints
{
	// start point
    double aCosEta1 = _a * cos(_eta1);
    double bSinEta1 = _b * sin(_eta1);
    _x1 = _cx + aCosEta1 * _cosTheta - bSinEta1 * _sinTheta;
    _y1 = _cy + aCosEta1 * _sinTheta + bSinEta1 * _cosTheta;
	
    // end point
    double aCosEta2 = _a * cos(_eta2);
    double bSinEta2 = _b * sin(_eta2);
    _x2 = _cx + aCosEta2 * _cosTheta - bSinEta2 * _sinTheta;
    _y2 = _cy + aCosEta2 * _sinTheta + bSinEta2 * _cosTheta;
}

-(void)computeBounds
{
	double bOnA = _b / _a;
    double etaXMin, etaXMax, etaYMin, etaYMax;
    if (abs(_sinTheta) < 0.1) {
		double tanTheta = _sinTheta / _cosTheta;
		if (_cosTheta < 0) {
			etaXMin = -atan(tanTheta * bOnA);
			etaXMax = etaXMin + M_PI;
			etaYMin = 0.5 * M_PI - atan(tanTheta / bOnA);
			etaYMax = etaYMin + M_PI;
		} else {
			etaXMax = -atan(tanTheta * bOnA);
			etaXMin = etaXMax - M_PI;
			etaYMax = 0.5 * M_PI - atan(tanTheta / bOnA);
			etaYMin = etaYMax - M_PI;
		}
    } else {
		double invTanTheta = _cosTheta / _sinTheta;
		if (_sinTheta < 0) {
			etaXMax = 0.5 * M_PI + atan(invTanTheta / bOnA);
			etaXMin = etaXMax - M_PI;
			etaYMin = atan(invTanTheta * bOnA);
			etaYMax = etaYMin + M_PI;
		} else {
			etaXMin = 0.5 * M_PI + atan(invTanTheta / bOnA);
			etaXMax = etaXMin + M_PI;
			etaYMax = atan(invTanTheta * bOnA);
			etaYMin = etaYMax - M_PI;
		}
    }
	
    etaXMin -= (2*M_PI) * floor((etaXMin - _eta1) / (2*M_PI));
    etaYMin -= (2*M_PI) * floor((etaYMin - _eta1) / (2*M_PI));
    etaXMax -= (2*M_PI) * floor((etaXMax - _eta1) / (2*M_PI));
    etaYMax -= (2*M_PI) * floor((etaYMax - _eta1) / (2*M_PI));
	
    _xLeft = (etaXMin <= _eta2)
	? (_cx + _a * cos(etaXMin) * _cosTheta - _b * sin(etaXMin) * _sinTheta)
	: MIN(_x1, _x2);
    _yUp = (etaYMin <= _eta2)
	? (_cy + _a * cos(etaYMin) * _sinTheta + _b * sin(etaYMin) * _cosTheta)
	: MIN(_y1, _y2);
    _width = ((etaXMax <= _eta2)
			  ? (_cx + _a * cos(etaXMax) * _cosTheta - _b * sin(etaXMax) * _sinTheta)
			  : MAX(_x1, _x2)) - _xLeft;
    _height = ((etaYMax <= _eta2)
			   ? (_cy + _a * cos(etaYMax) * _sinTheta + _b * sin(etaYMax) * _cosTheta)
			   : MAX(_y1, _y2)) - _yUp;
}

-(void)computeDerivedFlatnessParameters
{
	_f   = (_a - _b) / _a;
    _e2  = _f * (2.0 - _f);
    _g   = 1.0 - _f;
    _g2  = _g * _g;
}


#pragma mark - helper methods

+(double)rationalFunctionWithAbscissa:(double)x coefficients:(double [])c
{
	return (x * (x * c[0] + c[1]) + c[2]) / (x + c[3]);
}

-(double)estimateErrorWithDegree:(int)degree
					  startAngle:(double)etaA
						endAngle:(double)etaB
{
	double eta  = 0.5 * (etaA + etaB);
	
    if (degree < 2)
	{
		// start point
		double aCosEtaA  = _a * cos(etaA);
		double bSinEtaA  = _b * sin(etaA);
		double xA        = _cx + aCosEtaA * _cosTheta - bSinEtaA * _sinTheta;
		double yA        = _cy + aCosEtaA * _sinTheta + bSinEtaA * _cosTheta;
		
		// end point
		double aCosEtaB  = _a * cos(etaB);
		double bSinEtaB  = _b * sin(etaB);
		double xB        = _cx + aCosEtaB * _cosTheta - bSinEtaB * _sinTheta;
		double yB        = _cy + aCosEtaB * _sinTheta + bSinEtaB * _cosTheta;
		
		// maximal error point
		double aCosEta   = _a * cos(eta);
		double bSinEta   = _b * sin(eta);
		double x         = _cx + aCosEta * _cosTheta - bSinEta * _sinTheta;
		double y         = _cy + aCosEta * _sinTheta + bSinEta * _cosTheta;
		
		double dx = xB - xA;
		double dy = yB - yA;
		
		return abs(x * dy - y * dx + xB * yA - xA * yB) / sqrt(dx * dx + dy * dy);
    }
	else
	{
		double x    = _b / _a;
		double dEta = etaB - etaA;
		double cos2 = cos(2 * eta);
		double cos4 = cos(4 * eta);
		double cos6 = cos(6 * eta);
		
		// select the right coeficients set according to degree and _b/_a
		double (*coeffs)[4][4];
		double *safety;
		if (degree == 2) {
			coeffs = (x < 0.25) ? coeffs2Low : coeffs2High;
			safety = safety2;
		} else {
			coeffs = (x < 0.25) ? coeffs3Low : coeffs3High;
			safety = safety3;
		}
		
		double c0 = [EllipticalArc rationalFunctionWithAbscissa:x coefficients:coeffs[0][0]];
		c0 += cos2 * [EllipticalArc rationalFunctionWithAbscissa:x coefficients:coeffs[0][1]];
		c0 += cos4 * [EllipticalArc rationalFunctionWithAbscissa:x coefficients:coeffs[0][2]];
		c0 += cos6 * [EllipticalArc rationalFunctionWithAbscissa:x coefficients:coeffs[0][3]];
		
		double c1 = [EllipticalArc rationalFunctionWithAbscissa:x coefficients:coeffs[1][0]];
		c1 += cos2 * [EllipticalArc rationalFunctionWithAbscissa:x coefficients:coeffs[1][1]];
		c1 += cos4 * [EllipticalArc rationalFunctionWithAbscissa:x coefficients:coeffs[1][2]];
		c1 += cos6 * [EllipticalArc rationalFunctionWithAbscissa:x coefficients:coeffs[1][3]];
		
		return [EllipticalArc rationalFunctionWithAbscissa:x coefficients:safety] * _a * exp(c0 + c1 * dEta);
    }
}

-(CGPoint)pointAtAngle:(double)lambda
{
	CGPoint p;
	
	double eta      = atan2(sin(lambda) / _b, cos(lambda) / _a);
    double aCosEta  = _a * cos(eta);
    double bSinEta  = _b * sin(eta);
	
    p.x = _cx + aCosEta * _cosTheta - bSinEta * _sinTheta;
    p.y = _cy + aCosEta * _sinTheta + bSinEta * _cosTheta;
	
    return p;
}


#pragma mark - contains

-(BOOL)containsPoint:(CGPoint)point
{
	// position relative to the focii
    double dx1 = point.x - _xF1;
    double dy1 = point.y - _yF1;
    double dx2 = point.x - _xF2;
    double dy2 = point.y - _yF2;
    if ((dx1 * dx1 + dy1 * dy1 + dx2 * dx2 + dy2 * dy2) > (4 * _a * _a)) {
		// the point is outside of the ellipse
		return false;
    }
	
    if (_isPieSlice)
	{
		// check the location of the test point with respect to the
		// angular sector counted from the center of the ellipse
		double dxC = point.x - _cx;
		double dyC = point.y - _cy;
		double u   = dxC * _cosTheta + dyC * _sinTheta;
		double v   = dyC * _cosTheta - dxC * _sinTheta;
		double eta = atan2(v / _b, u / _a);
		eta -= (2*M_PI) * floor((eta - _eta1) / (2*M_PI));
		return (eta <= _eta2);
    }
	else {
		// check the location of the test point with respect to the
		// line joining the start and end points
		double dx = _x2 - _x1;
		double dy = _y2 - _y1;
		return ((point.x * dy - point.y * dx + _x2 * _y1 - _x1 * _y2) >= 0);
    }
}

-(BOOL)containsRectangle:(CGRect)rect
{
	double xPlusW = rect.origin.x + rect.size.width;
    double yPlusH = rect.origin.y + rect.size.height;
    return ([self containsPoint:rect.origin]
            && [self containsPoint:CGPointMake(xPlusW, rect.origin.y)]
            && [self containsPoint:CGPointMake(rect.origin.x, yPlusH)]
            && [self containsPoint:CGPointMake(xPlusW, yPlusH)]
            && (![self outlineIntersectLineWithStartPoint:rect.origin endPoint:CGPointMake(xPlusW, rect.origin.y)])
            && (![self outlineIntersectLineWithStartPoint:CGPointMake(xPlusW, rect.origin.y) endPoint:CGPointMake(xPlusW, yPlusH)])
            && (![self outlineIntersectLineWithStartPoint:CGPointMake(xPlusW, yPlusH) endPoint:CGPointMake(rect.origin.x, yPlusH)])
            && (![self outlineIntersectLineWithStartPoint:CGPointMake(rect.origin.x, yPlusH) endPoint:rect.origin]));
}

#pragma mark - intersects

-(BOOL)intersectLineWithStartPoint:(CGPoint)start endPoint:(CGPoint)end
{
	double dx = start.x - end.x;
    double dy = start.y - end.y;
    double l  = sqrt(dx * dx + dy * dy);
    if (l < (1.0e-10 * _a)) {
		// too small line segment, we consider it doesn't intersect anything
		return false;
    }
    double cz = (dx * _cosTheta + dy * _sinTheta) / l;
    double sz = (dy * _cosTheta - dx * _sinTheta) / l;
	
    // express position of the first point in canonical frame
    dx = start.x - _cx;
    dy = start.y - _cy;
    double u = dx * _cosTheta + dy * _sinTheta;
    double v = dy * _cosTheta - dx * _sinTheta;
	
    double u2         = u * u;
    double v2         = v * v;
    double g2u2ma2    = _g2 * (u2 - _a * _a);
	//    double g2u2ma2mv2 = g2u2ma2 - v2;
    double g2u2ma2pv2 = g2u2ma2 + v2;
	
    // compute intersections with the ellipse along the line
    // as the roots of _a 2nd degree polynom : c0 k^2 - 2 c1 k + c2 = 0
    double c0   = 1.0 - _e2 * cz * cz;
    double c1   = _g2 * u * cz + v * sz;
    double c2   = g2u2ma2pv2;
    double c12  = c1 * c1;
    double c0c2 = c0 * c2;
	
    if (c12 < c0c2) {
		// the line does not intersect the ellipse at all
		return false;
    }
	
    double k = (c1 >= 0)
	? (c1 + sqrt(c12 - c0c2)) / c0
	: c2 / (c1 - sqrt(c12 - c0c2));
    if ((k >= 0) && (k <= l)) {
		double uIntersect = u - k * cz;
		double vIntersect = v - k * sz;
		double eta = atan2(vIntersect / _b, uIntersect / _a);
		eta -= (2*M_PI) * floor((eta - _eta1) / (2*M_PI));
		if (eta <= _eta2) {
			return true;
		}
    }
	
    k = c2 / (k * c0);
    if ((k >= 0) && (k <= l)) {
		double uIntersect = u - k * cz;
		double vIntersect = v - k * sz;
		double eta = atan2(vIntersect / _b, uIntersect / _a);
		eta -= (2*M_PI) * floor((eta - _eta1) / (2*M_PI));
		if (eta <= _eta2) {
			return true;
		}
    }
	
    return false;
}

+(BOOL)lineIntersectLineWithStartPoint1:(CGPoint)start1
							  endPoint1:(CGPoint)end1
							startPoint2:(CGPoint)start2
							  endPoint2:(CGPoint)end2
{
	// elements of the equation of the (1, 2) line segment
    double dx12 = end1.x - start1.x;
    double dy12 = end1.y - start1.y;
    double k12  = end1.x * start1.y - start1.x * end1.y;
	
    // elements of the equation of the (A, B) line segment
    double dxAB = end2.x - start2.x;
    double dyAB = end2.y - start2.y;
    double kAB  = end2.x * start2.y - start2.x * end2.y;
	
    // compute relative positions of endpoints versus line segments
    double pAvs12 = start2.x * dy12 - start2.y * dx12 + k12;
    double pBvs12 = end2.x * dy12 - end2.y * dx12 + k12;
    double p1vsAB = start1.x * dyAB - start1.y * dxAB + kAB;
    double p2vsAB = start2.x * dyAB - start2.y * dxAB + kAB;
	
    return (pAvs12 * pBvs12 <= 0) && (p1vsAB * p2vsAB <= 0);
}

-(BOOL)outlineIntersectLineWithStartPoint:(CGPoint)start endPoint:(CGPoint)end
{
	if ([self intersectLineWithStartPoint:start endPoint:end])
		return true;
	
    if (_isPieSlice)
	{
		return [EllipticalArc lineIntersectLineWithStartPoint1:CGPointMake(_cx, _cy)
													 endPoint1:CGPointMake(_x1, _y1)
												   startPoint2:start endPoint2:end]
			|| [EllipticalArc lineIntersectLineWithStartPoint1:CGPointMake(_cx, _cy)
													 endPoint1:CGPointMake(_x1, _y1)
												   startPoint2:start endPoint2:end];
    }
	else {
		return [EllipticalArc lineIntersectLineWithStartPoint1:CGPointMake(_x1, _y1)
													 endPoint1:CGPointMake(_x2, _y2)
												   startPoint2:start endPoint2:end];
    }
}

-(BOOL)intersectsRectangle:(CGRect)rect
{
	double xPlusW = rect.origin.x + rect.size.width;
    double yPlusH = rect.origin.y + rect.size.height;
    return ([self containsPoint:rect.origin]
            || [self containsPoint:CGPointMake(xPlusW, rect.origin.y)]
            || [self containsPoint:CGPointMake(rect.origin.x, yPlusH)]
            || [self containsPoint:CGPointMake(xPlusW, yPlusH)]
            || (![self outlineIntersectLineWithStartPoint:rect.origin endPoint:CGPointMake(xPlusW, rect.origin.y)])
            || (![self outlineIntersectLineWithStartPoint:CGPointMake(xPlusW, rect.origin.y) endPoint:CGPointMake(xPlusW, yPlusH)])
            || (![self outlineIntersectLineWithStartPoint:CGPointMake(xPlusW, yPlusH) endPoint:CGPointMake(rect.origin.x, yPlusH)])
            || (![self outlineIntersectLineWithStartPoint:CGPointMake(rect.origin.x, yPlusH) endPoint:rect.origin]));
}


#pragma mark - builders

-(CGRect)bounds
{
	return CGRectMake(_xLeft, _yUp, _width, _height);
}

-(CGMutablePathRef)pathWithDegree:(int)degree
						threshold:(double)threshold
					tranformation:(CGAffineTransform *)m
{
	// find the number of Bézier curves needed
    BOOL found = NO;
    int n = 1;
    while ((!found) && (n < 1024))
	{
		double dEta = (_eta2 - _eta1) / n;
		if (dEta <= 0.5 * M_PI)
		{
			double etaB = _eta1;
			found = YES;
			for (int i = 0; found && (i < n); ++i)
			{
				double etaA = etaB;
				etaB += dEta;
				found = ([self estimateErrorWithDegree:degree startAngle:etaA endAngle:etaB] <= threshold);
			}
		}
		n = n << 1;
    }
	
    CGMutablePathRef path = CGPathCreateMutable();
    double dEta = (_eta2 - _eta1) / n;
    double etaB = _eta1;
	
    double cosEtaB  = cos(etaB);
    double sinEtaB  = sin(etaB);
    double aCosEtaB = _a * cosEtaB;
    double bSinEtaB = _b * sinEtaB;
    double aSinEtaB = _a * sinEtaB;
    double bCosEtaB = _b * cosEtaB;
    double xB       = _cx + aCosEtaB * _cosTheta - bSinEtaB * _sinTheta;
    double yB       = _cy + aCosEtaB * _sinTheta + bSinEtaB * _cosTheta;
    double xBDot    = -aSinEtaB * _cosTheta - bCosEtaB * _sinTheta;
    double yBDot    = -aSinEtaB * _sinTheta + bCosEtaB * _cosTheta;
	
    if (_isPieSlice)
	{
		CGPathMoveToPoint(path, NULL, _cx, _cy);
		CGPathAddLineToPoint(path, NULL, xB, yB);
    }
	else CGPathMoveToPoint(path, NULL, xB, yB);
	
    double t     = tan(0.5 * dEta);
    double alpha = sin(dEta) * (sqrt(4 + 3 * t * t) - 1) / 3;
	
    for (int i = 0; i < n; ++i)
	{
//		double etaA  = etaB;
		double xA    = xB;
		double yA    = yB;
		double xADot = xBDot;
		double yADot = yBDot;
		
		etaB    += dEta;
		cosEtaB  = cos(etaB);
		sinEtaB  = sin(etaB);
		aCosEtaB = _a * cosEtaB;
		bSinEtaB = _b * sinEtaB;
		aSinEtaB = _a * sinEtaB;
		bCosEtaB = _b * cosEtaB;
		xB       = _cx + aCosEtaB * _cosTheta - bSinEtaB * _sinTheta;
		yB       = _cy + aCosEtaB * _sinTheta + bSinEtaB * _cosTheta;
		xBDot    = -aSinEtaB * _cosTheta - bCosEtaB * _sinTheta;
		yBDot    = -aSinEtaB * _sinTheta + bCosEtaB * _cosTheta;
		
		if (degree == 1)
		{
			CGPathMoveToPoint(path, NULL, xB, yB);
		}
		else if (degree == 2)
		{
			double k = (yBDot * (xB - xA) - xBDot * (yB - yA)) / (xADot * yBDot - yADot * xBDot);
			CGPathAddQuadCurveToPoint(path, NULL, (xA + k * xADot), (yA + k * yADot), xB, yB);
		}
		else {
			CGPathAddCurveToPoint(path, NULL, (xA + alpha * xADot), (yA + alpha * yADot), (xB - alpha * xBDot), (yB - alpha * yBDot), xB, yB);
		}
    }
	
    if (_isPieSlice)
	{
		CGPathCloseSubpath(path);
    }
	
	CGMutablePathRef returnPath = CGPathCreateMutableCopyByTransformingPath(path, m);
	CGPathRelease(path);
	
	return returnPath;
}

-(CGMutablePathRef)pathWithTranformation:(CGAffineTransform *)m
{
	return [self pathWithDegree:_maxDegree threshold:_defaultFlatness tranformation:m];
}

-(CGMutablePathRef)pathWithTranformation:(CGAffineTransform *)m flatness:(double)flatness
{
	return [self pathWithDegree:_maxDegree threshold:flatness tranformation:m];
}

@end
