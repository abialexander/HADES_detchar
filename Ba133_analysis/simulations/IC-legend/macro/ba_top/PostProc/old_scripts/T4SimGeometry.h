
//detector geometry

#ifndef GADA_T4SimGeometry_H
#define GADA_T4SimGeometry_H

//root
#include "TF1.h"

//gerda-ada
//#include "T4SimHitAnalysis.h"
//#include "T4SimTransitionLayer.h"
//c++
#include <algorithm>

//namespace gada {

struct T4TwoVector
{
   double x,y;

   T4TwoVector& operator+= (const T4TwoVector& rhs) {
      x += rhs.x; y += rhs.y;
      return *this;
   };
   T4TwoVector& operator-= (const T4TwoVector& rhs) {
      x -= rhs.x; y -= rhs.y;
      return *this;
   };
   T4TwoVector  operator+  (const T4TwoVector& v) { return {x+v.x,y+v.y}; };
   T4TwoVector  operator-  (const T4TwoVector& v) { return {x-v.x,y-v.y}; };
   double       operator*  (const T4TwoVector& v) { return x*v.x+y*v.y;   };
   T4TwoVector  operator*  (const double rhs)     { return {x*rhs,y*rhs}; };
   T4TwoVector  operator/  (const double rhs)     { return {x/rhs,y/rhs}; };

   double length2() { return *this * *this;   };
   double length()  { return sqrt(length2()); };

//   friend std::ostream& operator<<(std::ostream& os, const T4TwoVector& v) {
 //     os << "(" << v.x << "," << v.y << ")";
 //     return os;
//   };
};

struct T4TwoDLine
{
   T4TwoVector p1,p2;

   double length()  { return (this->p1-this->p2).length();  };
   double length2() { return (this->p1-this->p2).length2(); };
   double distance(T4TwoVector v) {
      T4TwoVector proj = p1+(p2-p1)*std::max(0.,std::min((v-p1)*(p2-p1)/length2(),1.));
      return (proj-v).length();
   };

//   friend std::ostream& operator<<(std::ostream& os, const T4TwoDLine& l) {
 //     os << l.p1 << "--" << l.p2;
//      return os;
 //  };
};

class T4SimHPGe //(tappered) zylinder with groove (and borehole)
{
   private:
      std::string             fName;
      double                  fRadius;
      double                  fHeight;
      double                  fFCCD;
      double                  fDLT;
      std::vector<T4TwoDLine> fNPlus;
      std::vector<T4TwoDLine> fBore;
//      std::vector<T4TwoDLine> fPPlus;
//      std::vector<T4TwoDLine> fGroove;

   public:
      T4SimHPGe(std::string name,
       double fccd=0 , double dlt=0,
       double radius=0, double height=0,
       double grooveOuterRadius = 0,
       double grooveInnerRadius = 0,
       double grooveDepth = 0,
       double coneRadius  = 0,
       double coneHeight  = 0,
       double boreRadius  = 0,
       double boreDepth   = 0) {
        fName = name;
        fRadius = radius;
        fHeight = height;
	fFCCD= fccd;
	fDLT=dlt;
        if(!coneHeight) {
            fNPlus = {
              {{grooveOuterRadius,0.},{radius,0.}},  //bottom
              {{radius,0.},{radius,height}},         //side
              {{radius,height},{boreRadius,height}}, //top
            };
	    fBore={  
	      {{boreRadius,height},{boreRadius,height-boreDepth}},  //top bore hole
              {{boreRadius,height-boreDepth},{0.,height-boreDepth}} //top bore hole
            };
	}
        else {//tappered
            fNPlus = {
              {{grooveOuterRadius,0.},{radius,0.}},             //bottom
              {{radius,0.},{radius,height-coneHeight}},         //side
              {{radius,height-coneHeight},{coneRadius,height}}, //tapper
              {{coneRadius,height},{boreRadius,height}},        //top
            };
	    fBore={  
	      {{boreRadius,height},{boreRadius,height-boreDepth}},  //top bore hole
              {{boreRadius,height-boreDepth},{0.,height-boreDepth}} //top bore hole
            };
	}
/*        fPPlus = {
            {{grooveInnerRadius,0.},{0.,0.}} //outer
          };
        fGroove = {
          {{grooveInnerRadius,0.},{grooveInnerRadius,grooveDepth}},
          {{grooveInnerRadius,grooveDepth},{grooveOuterRadius,grooveDepth}},
          {{grooveOuterRadius,grooveDepth},{grooveOuterRadius,0.}}
        };
*/

     };
      virtual ~T4SimHPGe() {};

      std::string              GetName()        { return fName;         }
      double                   GetRadius()      { return  fRadius;      }
      double                   GetHeight()      { return  fHeight;      }
      double                   GetFCCD()        { return  fFCCD;        }
      double                   GetDLT()         { return  fDLT;         }
      std::vector<T4TwoDLine>* GetNPlus()       { return &fNPlus;       }
      std::vector<T4TwoDLine>* GetBore()        { return &fBore;        }
//      std::vector<T4TwoDLine>* GetPPlus()       { return &fPPlus;       }
//      std::vector<T4TwoDLine>* GetGroove()      { return &fGroove;      }



      void SetName(std::string name)                        { fName   = name;           }
      void SetHeight(double height)                         { fHeight=height;           }
      void SetRadius(double radius)                         { fRadius=radius;           }
      void SetFCCD(double fccd)                             { fFCCD=fccd;               }
      void SetDLT(double dlt)                               { fDLT=dlt;                 }
      void SetNPlus(std::vector<T4TwoDLine> nPlus)          { fNPlus=nPlus;             }
      void SetBore(std::vector<T4TwoDLine> bore)            { fBore=bore;               }
//      void SetPPlus(std::vector<T4TwoDLine> pPlus)          { fPPlus=pPlus;             }
//      void SetGroove(std::vector<T4TwoDLine> groove)        { fGroove=groove;           }

	
	double  FCCDBore(double x){
		if 	  	(x <= fDLT/2) 			   	     return 0.*x;
    		else if 	(fDLT!=fFCCD && x>fDLT/2. && x<fFCCD/2.)     return 2./(fFCCD-fDLT)*x-fDLT/(fFCCD-fDLT);
    		else						   	     return 1.+0.*x; 
	}

	double  FCCDOuter(double x){
		if 	  	(x <= fDLT) 			  	return 0.*x;
    		else if 	(fDLT!=fFCCD && x>fDLT && x<fFCCD)      return 1./(fFCCD-fDLT)*x-fDLT/(fFCCD-fDLT);
    		else						  	return 1.+0.*x;
	}


      double GetChargeCollectionEfficiency(double radius,double z) {
        double distanceToNPlus  = GetDistanceToNPlus(radius,z);
        double distanceToBore   = GetDistanceToBore(radius,z);
//        double distanceToPPlus  = GetDistanceToPPlus(radius,z);
//        double distanceToGroove = GetDistanceToGroove(radius,z);

        // check hit inside detector
	 //double minDist = std::min( std::min ( std::min(distanceToBore, distanceToNPlus), distanceToPPlus), distanceToGroove);
	 double minDist = std::min( distanceToBore, distanceToNPlus);
        
	if( (radius>fRadius || z>fHeight || z<0) && minDist > 2e-5 ) {
          std::cerr << "error : in " << fName << " hit outside detector" << std::endl;
          std::cerr << "\tmin distance to contact : " << minDist << std::endl;
          std::cerr << "\tz : " << z << std::endl;
          std::cerr << "\tr : " << radius << std::endl;
        }

        if(minDist < 0) return 0;
	

	else if(minDist==distanceToBore){
		return FCCDBore(minDist) ;
	}
	else{ 
		return FCCDOuter(minDist);
	};



      }

   protected:
      double GetMinimumDistance(std::vector<T4TwoDLine> chain,T4TwoVector point) {
        if(!chain.size()) return 0.;
          double distance = chain.at(0).distance(point);
          for(int i=1;i<(int)chain.size();i++)
          distance =  std::min(distance,chain.at(i).distance(point));
        return distance;
      }
      double GetDistanceToNPlus(double radius,double z) { //only n+
        return GetMinimumDistance(fNPlus,T4TwoVector({radius,z}));
      }
      double GetDistanceToBore(double radius,double z) { //only borehole
        return GetMinimumDistance(fBore,T4TwoVector({radius,z}));
      }
/*      double GetDistanceToPPlus(double radius,double z) { //only p+
        return GetMinimumDistance(fPPlus,T4TwoVector({radius,z}));
      }
      double GetDistanceToGroove(double radius,double z) { //only groove
        return GetMinimumDistance(fGroove,T4TwoVector({radius,z}));
      }
*/
};
/*
class T4SimPlacedHPGe : public T4SimHPGe
{
   public:
      T4SimPlacedHPGe(T4SimHPGe hpge,T4ThreeVector position,bool bottomUp=false) :
       T4SimHPGe(hpge), fPosition(position), fBottomUp(bottomUp) {
      };
      virtual ~T4SimPlacedHPGe() {};

      T4ThreeVector GetPosition() { return  fPosition; };
      bool          GetBottomUp() { return  fBottomUp; };

      void FlipMe(bool bottomUp=true) { fBottomUp=bottomUp; };

      double GetChargeCollectionEfficiency(T4ThreeVector hitPosition) {
        hitPosition-=fPosition; //rel position
        return T4SimHPGe::GetChargeCollectionEfficiency(hitPosition.radius(),
          fBottomUp ? GetHeight()-hitPosition.z : hitPosition.z);
      };

   private:
      T4ThreeVector fPosition; //p+ center position
      bool          fBottomUp; //face p+ up
};
*/
//} // namespace gada

#endif
