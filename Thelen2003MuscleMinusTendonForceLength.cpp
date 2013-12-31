// Thelen2003MuscleMinusTendonForceLength.cpp
/*
 * Copyright (c)  2006, Stanford University. All rights reserved. 
* Use of the OpenSim software in source form is permitted provided that the following
* conditions are met:
* 	1. The software is used only for non-commercial research and education. It may not
*     be used in relation to any commercial activity.
* 	2. The software is not distributed or redistributed.  Software distribution is allowed 
*     only through https://simtk.org/home/opensim.
* 	3. Use of the OpenSim software or derivatives must be acknowledged in all publications,
*      presentations, or documents describing work in which OpenSim or derivatives are used.
* 	4. Credits to developers may not be removed from executables
*     created from modifications of the source.
* 	5. Modifications of source code must retain the above copyright notice, this list of
*     conditions and the following disclaimer. 
* 
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
*  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
*  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
*  SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
*  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
*  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
*  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
*  OR BUSINESS INTERRUPTION) OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
*  WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

//=============================================================================
// INCLUDES
//=============================================================================
#include "Thelen2003MuscleMinusTendonForceLength.h"
#include <OpenSim/Common/SimmMacros.h>
#include <OpenSim/Common/rdMath.h>
#include <OpenSim/Common/DebugUtilities.h>

//=============================================================================
// STATICS
//=============================================================================
using namespace std;
using namespace OpenSim;

const int Thelen2003MuscleMinusTendonForceLength::STATE_ACTIVATION = 0;
const int Thelen2003MuscleMinusTendonForceLength::STATE_FIBER_LENGTH = 1;

//=============================================================================
// CONSTRUCTOR(S) AND DESTRUCTOR
//=============================================================================
//_____________________________________________________________________________
/**
 * Default constructor.
 */
Thelen2003MuscleMinusTendonForceLength::Thelen2003MuscleMinusTendonForceLength() :
   AbstractMuscle(),
	_maxIsometricForce(_maxIsometricForceProp.getValueDbl()),
	_optimalFiberLength(_optimalFiberLengthProp.getValueDbl()),
	_tendonSlackLength(_tendonSlackLengthProp.getValueDbl()),
	_pennationAngle(_pennationAngleProp.getValueDbl()),
	_activationTimeConstant(_activationTimeConstantProp.getValueDbl()),
	_deactivationTimeConstant(_deactivationTimeConstantProp.getValueDbl()),
	_vmax(_vmaxProp.getValueDbl()),
	_vmax0(_vmax0Prop.getValueDbl()),
	_fmaxTendonStrain(_fmaxTendonStrainProp.getValueDbl()),
	_fmaxMuscleStrain(_fmaxMuscleStrainProp.getValueDbl()),
	_kShapeActive(_kShapeActiveProp.getValueDbl()),
	_kShapePassive(_kShapePassiveProp.getValueDbl()),
	_damping(_dampingProp.getValueDbl()),
	_af(_afProp.getValueDbl()),
	_flen(_flenProp.getValueDbl())
{
	setNull();
	setupProperties();
}

//_____________________________________________________________________________
/**
 * Destructor.
 */
Thelen2003MuscleMinusTendonForceLength::~Thelen2003MuscleMinusTendonForceLength()
{
}

//_____________________________________________________________________________
/**
 * Copy constructor.
 *
 * @param aMuscle Thelen2003MuscleMinusTendonForceLength to be copied.
 */
Thelen2003MuscleMinusTendonForceLength::Thelen2003MuscleMinusTendonForceLength(const Thelen2003MuscleMinusTendonForceLength &aMuscle) :
   AbstractMuscle(aMuscle),
	_maxIsometricForce(_maxIsometricForceProp.getValueDbl()),
	_optimalFiberLength(_optimalFiberLengthProp.getValueDbl()),
	_tendonSlackLength(_tendonSlackLengthProp.getValueDbl()),
	_pennationAngle(_pennationAngleProp.getValueDbl()),
	_activationTimeConstant(_activationTimeConstantProp.getValueDbl()),
	_deactivationTimeConstant(_deactivationTimeConstantProp.getValueDbl()),
	_vmax(_vmaxProp.getValueDbl()),
	_vmax0(_vmax0Prop.getValueDbl()),
	_fmaxTendonStrain(_fmaxTendonStrainProp.getValueDbl()),
	_fmaxMuscleStrain(_fmaxMuscleStrainProp.getValueDbl()),
	_kShapeActive(_kShapeActiveProp.getValueDbl()),
	_kShapePassive(_kShapePassiveProp.getValueDbl()),
	_damping(_dampingProp.getValueDbl()),
	_af(_afProp.getValueDbl()),
	_flen(_flenProp.getValueDbl())
{
	setNull();
	setupProperties();
	copyData(aMuscle);
	setup(aMuscle.getModel());
}

//_____________________________________________________________________________
/**
 * Copy this muscle point and return a pointer to the copy.
 * The copy constructor for this class is used.
 *
 * @return Pointer to a copy of this Thelen2003MuscleMinusTendonForceLength.
 */
Object* Thelen2003MuscleMinusTendonForceLength::copy() const
{
	Thelen2003MuscleMinusTendonForceLength *musc = new Thelen2003MuscleMinusTendonForceLength(*this);
	return(musc);
}

//=============================================================================
// CONSTRUCTION METHODS
//=============================================================================
//_____________________________________________________________________________
/**
 * Copy data members from one Thelen2003MuscleMinusTendonForceLength to another.
 *
 * @param aMuscle Thelen2003MuscleMinusTendonForceLength to be copied.
 */
void Thelen2003MuscleMinusTendonForceLength::copyData(const Thelen2003MuscleMinusTendonForceLength &aMuscle)
{
	_maxIsometricForce = aMuscle._maxIsometricForce;
	_optimalFiberLength = aMuscle._optimalFiberLength;
	_tendonSlackLength = aMuscle._tendonSlackLength;
	_pennationAngle = aMuscle._pennationAngle;
	_activationTimeConstant = aMuscle._activationTimeConstant;
	_deactivationTimeConstant = aMuscle._deactivationTimeConstant;
	_vmax = aMuscle._vmax;
	_vmax0 = aMuscle._vmax0;
	_fmaxTendonStrain = aMuscle._fmaxTendonStrain;
	_fmaxMuscleStrain = aMuscle._fmaxMuscleStrain;
	_kShapeActive = aMuscle._kShapeActive;
	_kShapePassive = aMuscle._kShapePassive;
	_damping = aMuscle._damping;
	_af = aMuscle._af;
	_flen = aMuscle._flen;
}

//_____________________________________________________________________________
/**
 * Set the data members of this Thelen2003MuscleMinusTendonForceLength to their null values.
 */
void Thelen2003MuscleMinusTendonForceLength::setNull()
{
	setType("Thelen2003MuscleMinusTendonForceLength");

	setNumControls(1); setNumStates(2); setNumPseudoStates(0);
	bindControl(0, _excitation, "excitation");
	bindState(STATE_ACTIVATION, _activation, "activation");
	bindState(STATE_FIBER_LENGTH, _fiberLength, "fiber_length");

	_excitation = 0.0;
	_activation = 0.0;
	_fiberLength = 0.0;
}

//_____________________________________________________________________________
/**
 * Connect properties to local pointers.
 */
void Thelen2003MuscleMinusTendonForceLength::setupProperties()
{
	_maxIsometricForceProp.setName("max_isometric_force");
	_maxIsometricForceProp.setValue(546.0);
   _maxIsometricForceProp.setComment("maximum isometric force of the muscle fibers");
	_propertySet.append(&_maxIsometricForceProp, "Parameters");

	_optimalFiberLengthProp.setName("optimal_fiber_length");
	_optimalFiberLengthProp.setValue(0.0535);
	_optimalFiberLengthProp.setComment("optimal length of the muscle fibers");
	_propertySet.append(&_optimalFiberLengthProp, "Parameters");

	_tendonSlackLengthProp.setName("tendon_slack_length");
	_tendonSlackLengthProp.setValue(0.078);
	_tendonSlackLengthProp.setComment("resting length of the tendon");
	_propertySet.append(&_tendonSlackLengthProp, "Parameters");

	_pennationAngleProp.setName("pennation_angle");
	_pennationAngleProp.setComment("angle between tendon and fibers at optimal fiber length");
	_pennationAngleProp.setValue(0.0);
	_propertySet.append(&_pennationAngleProp, "Parameters");

	_activationTimeConstantProp.setName("activation_time_constant");
	_activationTimeConstantProp.setValue(0.01);
	_activationTimeConstantProp.setComment("time constant for ramping up of muscle activation");
	_propertySet.append(&_activationTimeConstantProp, "Parameters");

	_deactivationTimeConstantProp.setName("deactivation_time_constant");
	_deactivationTimeConstantProp.setValue(0.04);
	_deactivationTimeConstantProp.setComment("time constant for ramping down of muscle activation");
	_propertySet.append(&_deactivationTimeConstantProp, "Parameters");

	_vmaxProp.setName("Vmax");
	_vmaxProp.setValue(10.0);
	_vmaxProp.setComment("maximum contraction velocity at full activation in fiber lengths per second");
	_propertySet.append(&_vmaxProp, "Parameters");

	_vmax0Prop.setName("Vmax0");
	_vmax0Prop.setValue(5.0);
	_vmax0Prop.setComment("maximum contraction velocity at low activation in fiber lengths per second");
	_propertySet.append(&_vmax0Prop, "Parameters");

	_fmaxTendonStrainProp.setName("FmaxTendonStrain");
	_fmaxTendonStrainProp.setValue(0.033);
	_fmaxTendonStrainProp.setComment("tendon strain due to maximum isometric muscle force");
	_propertySet.append(&_fmaxTendonStrainProp, "Parameters");

	_fmaxMuscleStrainProp.setName("FmaxMuscleStrain");
	_fmaxMuscleStrainProp.setValue(0.6);
	_fmaxMuscleStrainProp.setComment("passive muscle strain due to maximum isometric muscle force");
	_propertySet.append(&_fmaxMuscleStrainProp, "Parameters");

	_kShapeActiveProp.setName("KshapeActive");
	_kShapeActiveProp.setValue(0.5);
	_kShapeActiveProp.setComment("shape factor for Gaussian active muscle force-length relationship");
	_propertySet.append(&_kShapeActiveProp, "Parameters");

	_kShapePassiveProp.setName("KshapePassive");
	_kShapePassiveProp.setValue(4.0);
	_kShapePassiveProp.setComment("exponential shape factor for passive force-length relationship");
	_propertySet.append(&_kShapePassiveProp, "Parameters");

	_dampingProp.setName("damping");
	_dampingProp.setValue(0.05);
	_dampingProp.setComment("passive damping in the force-velocity relationship");
	_propertySet.append(&_dampingProp, "Parameters");

	_afProp.setName("Af");
	_afProp.setValue(0.3);
	_afProp.setComment("force-velocity shape factor");
	_propertySet.append(&_afProp, "Parameters");

	_flenProp.setName("Flen");
	_flenProp.setValue(1.8);
	_flenProp.setComment("maximum normalized lengthening force");
	_propertySet.append(&_flenProp, "Parameters");
}

//_____________________________________________________________________________
/**
 * Perform some set up functions that happen after the
 * object has been deserialized or copied.
 *
 * @param aModel model containing this Thelen2003MuscleMinusTendonForceLength.
 */
void Thelen2003MuscleMinusTendonForceLength::setup(Model* aModel)
{
	// Base class
	AbstractMuscle::setup(aModel);

	// aModel will be NULL when objects are being registered.
	if (aModel == NULL)
		return;

	// Reasonable initial activation value
	_activation = 0.01;

	// Compute isometric force to get starting value
	// of _fiberLength.
   computeEquilibrium();
}

//_____________________________________________________________________________
/**
 * Copy the property values from another actuator, which may not be
 * a Thelen2003MuscleMinusTendonForceLength.
 *
 * @param aActuator Actuator to copy property values from.
 */
void Thelen2003MuscleMinusTendonForceLength::copyPropertyValues(AbstractActuator& aActuator)
{
	AbstractMuscle::copyPropertyValues(aActuator);

	const Property* prop = aActuator.getPropertySet().contains("max_isometric_force");
	if (prop) _maxIsometricForceProp.setValue(prop->getValueDbl());

	prop = aActuator.getPropertySet().contains("optimal_fiber_length");
	if (prop) _optimalFiberLengthProp.setValue(prop->getValueDbl());

	prop = aActuator.getPropertySet().contains("tendon_slack_length");
	if (prop) _tendonSlackLengthProp.setValue(prop->getValueDbl());

	prop = aActuator.getPropertySet().contains("pennation_angle");
	if (prop) _pennationAngleProp.setValue(prop->getValueDbl());

	prop = aActuator.getPropertySet().contains("activation_time_constant");
	if (prop) _activationTimeConstantProp.setValue(prop->getValueDbl());

	prop = aActuator.getPropertySet().contains("deactivation_time_constant");
	if (prop) _deactivationTimeConstantProp.setValue(prop->getValueDbl());

	prop = aActuator.getPropertySet().contains("Vmax");
	if (prop) _vmaxProp.setValue(prop->getValueDbl());

	prop = aActuator.getPropertySet().contains("Vmax0");
	if (prop) _vmax0Prop.setValue(prop->getValueDbl());

	prop = aActuator.getPropertySet().contains("FmaxTendonStrain");
	if (prop) _fmaxTendonStrainProp.setValue(prop->getValueDbl());

	prop = aActuator.getPropertySet().contains("FmaxMuscleStrain");
	if (prop) _fmaxMuscleStrainProp.setValue(prop->getValueDbl());

	prop = aActuator.getPropertySet().contains("KshapeActive");
	if (prop) _kShapeActiveProp.setValue(prop->getValueDbl());

	prop = aActuator.getPropertySet().contains("KshapePassive");
	if (prop) _kShapePassiveProp.setValue(prop->getValueDbl());

	prop = aActuator.getPropertySet().contains("damping");
	if (prop) _dampingProp.setValue(prop->getValueDbl());

	prop = aActuator.getPropertySet().contains("Af");
	if (prop) _afProp.setValue(prop->getValueDbl());

	prop = aActuator.getPropertySet().contains("Flen");
	if (prop) _flenProp.setValue(prop->getValueDbl());
}

//=============================================================================
// OPERATORS
//=============================================================================
//_____________________________________________________________________________
/**
 * Assignment operator.
 *
 * @return Reference to this object.
 */
Thelen2003MuscleMinusTendonForceLength& Thelen2003MuscleMinusTendonForceLength::operator=(const Thelen2003MuscleMinusTendonForceLength &aMuscle)
{
	// BASE CLASS
	AbstractMuscle::operator=(aMuscle);

	copyData(aMuscle);

	setup(aMuscle.getModel());

	return(*this);
}


//=============================================================================
// SCALING
//=============================================================================
//_____________________________________________________________________________
/**
 * Scale the muscle.
 *
 * @param aScaleSet XYZ scale factors for the bodies
 * @return Whether or not the muscle was scaled successfully
 */
void Thelen2003MuscleMinusTendonForceLength::scale(const ScaleSet& aScaleSet)
{
	AbstractMuscle::scale(aScaleSet);

	// some force-generating parameters are scaled in postScale(),
	// so as of now there is nothing else to do here...
}

//_____________________________________________________________________________
/**
 * Perform computations that need to happen after the muscle is scaled.
 * For this object, that entails comparing the musculotendon length
 * before and after scaling, and scaling some of the force-generating
 * properties a proportional amount.
 *
 * @param aScaleSet XYZ scale factors for the bodies.
 */
void Thelen2003MuscleMinusTendonForceLength::postScale(const ScaleSet& aScaleSet)
{
	// Base class
	AbstractMuscle::postScale(aScaleSet);

	if (_preScaleLength > 0.0)
	{
		double scaleFactor = getLength() / _preScaleLength;

		_optimalFiberLength *= scaleFactor;
		_tendonSlackLength *= scaleFactor;
		//_maxIsometricForce *= scaleFactor;

		_preScaleLength = 0.0;
	}
}


//=============================================================================
// GET
//=============================================================================
//-----------------------------------------------------------------------------
// PENNATION ANGLE
//-----------------------------------------------------------------------------
//_____________________________________________________________________________
/**
 * Get the current pennation angle of the muscle fiber(s).
 *
 * @param Pennation angle.
 */
double Thelen2003MuscleMinusTendonForceLength::getPennationAngle()
{
	return calcPennation(_fiberLength,_optimalFiberLength,_pennationAngle);
}

//-----------------------------------------------------------------------------
// LENGTH
//-----------------------------------------------------------------------------
//_____________________________________________________________________________
/**
 * Get the length of the muscle fiber(s).
 *
 * @param Current length of the muscle fiber(s).
 */
double Thelen2003MuscleMinusTendonForceLength::getFiberLength()
{
	return _fiberLength;
}
//_____________________________________________________________________________
/**
 * Get the normalized length of the muscle fiber(s).  This is the current
 * fiber length(s) divided by the optimal fiber length.
 *
 * @param Current length of the muscle fiber(s).
 */
double Thelen2003MuscleMinusTendonForceLength::getNormalizedFiberLength()
{
	return _fiberLength / getOptimalFiberLength();
}

//-----------------------------------------------------------------------------
// FORCE
//-----------------------------------------------------------------------------
//_____________________________________________________________________________
/**
 * Get the passive force generated by the muscle fibers.
 *
 * @param Current passive force of the muscle fiber(s).
 */
double Thelen2003MuscleMinusTendonForceLength::getPassiveFiberForce()
{
	return _passiveForce;
}
//_____________________________________________________________________________
/**
 * Get the applied active force generated by the muscle fibers.
 * Thelen2003MuscleMinusTendonForceLength implements this function to apply the
 * active force that was actually computed instead of just subtracting the
 * passive fiber force from the total fiber force.
 *
 * @return Current active force of the muscle fibers.
 */
double Thelen2003MuscleMinusTendonForceLength::getAppliedActiveFiberForce()
{
	return _activeForce;
}


//=============================================================================
// COMPUTATION
//=============================================================================
//_____________________________________________________________________________
/**
 * Compute the derivatives of the muscle states.
 *
 * @param rDYDT the state derivatives are returned here.
 */
void Thelen2003MuscleMinusTendonForceLength::computeStateDerivatives(double rDYDT[])
{
	if (!rDYDT)
		return;

	rDYDT[STATE_ACTIVATION] = _activationDeriv;
	rDYDT[STATE_FIBER_LENGTH] = _fiberLengthDeriv;
}

//_____________________________________________________________________________
/**
 * Compute the equilibrium states.  This method computes a fiber length
 * for the muscle that is consistent with the muscle's activation level.
 */
void Thelen2003MuscleMinusTendonForceLength::computeEquilibrium()
{
	double force = computeIsometricForce(_activation);

	//cout<<getName()<<": isometric force = "<<force<<endl;
	//cout<<getName()<<": fiber length = "<<_fiberLength<<endl;
}

//_____________________________________________________________________________
/**
 * Compute the actuation for the muscle. This function assumes
 * that computeDerivatives has already been called.
 *
 * This function is based on muscle_deriv_func9 from derivs.c (old pipeline code)
 */
void Thelen2003MuscleMinusTendonForceLength::computeActuation()
{
	// Base Class (to calculate speed)
	AbstractMuscle::computeActuation();

   double normState[2], normStateDeriv[2], norm_tendon_length, ca;
   double norm_muscle_tendon_length, pennation_angle;

   /* Normalize the muscle states */
   normState[STATE_ACTIVATION] = _activation;
   normState[STATE_FIBER_LENGTH] = _fiberLength / _optimalFiberLength;

	// Maximum contraction velocity is an activation scaled value
	double Vmax = _vmax;
	if (normState[STATE_ACTIVATION]<1.0)
		Vmax = _vmax0 + normState[STATE_ACTIVATION]*(Vmax-_vmax0);
	Vmax = Vmax*_optimalFiberLength;

   /* Compute normalized muscle state derivatives */
   if (_excitation >= normState[STATE_ACTIVATION])
      normStateDeriv[STATE_ACTIVATION] = (_excitation - normState[STATE_ACTIVATION]) / _activationTimeConstant;
   else
      normStateDeriv[STATE_ACTIVATION] = (_excitation - normState[STATE_ACTIVATION]) / _deactivationTimeConstant;

   norm_muscle_tendon_length = getLength() / _optimalFiberLength;
   //norm_tendon_length = norm_muscle_tendon_length - normState[STATE_FIBER_LENGTH] * ca;
   norm_tendon_length = -1;
	string name = getName();
	if( name=="glut_med1_r" ) norm_tendon_length = 0.080208;
	else if( name=="glut_med2_r" ) norm_tendon_length = 0.054737;
	else if( name=="glut_med3_r" ) norm_tendon_length = 0.055189;
	else if( name=="glut_min1_r" ) norm_tendon_length = 0.016771;
	else if( name=="glut_min2_r" ) norm_tendon_length = 0.027448;
	else if( name=="glut_min3_r" ) norm_tendon_length = 0.053949;
	else if( name=="semimem_r" ) norm_tendon_length = 0.40048;
	else if( name=="semiten_r" ) norm_tendon_length = 0.28044;
	else if( name=="bifemlh_r" ) norm_tendon_length = 0.3614;
	else if( name=="bifemsh_r" ) norm_tendon_length = 0.097138;
	else if( name=="sar_r" ) norm_tendon_length = 0.10816;
	else if( name=="add_long_r" ) norm_tendon_length = 0.12198;
	else if( name=="add_brev_r" ) norm_tendon_length = 0.02184;
	else if( name=="add_mag1_r" ) norm_tendon_length = 0.066238;
	else if( name=="add_mag2_r" ) norm_tendon_length = 0.13609;
	else if( name=="add_mag3_r" ) norm_tendon_length = 0.28173;
	else if( name=="tfl_r" ) norm_tendon_length = 0.46172;
	else if( name=="pect_r" ) norm_tendon_length = 0.035827;
	else if( name=="grac_r" ) norm_tendon_length = 0.1386;
	else if( name=="glut_max1_r" ) norm_tendon_length = 0.12998;
	else if( name=="glut_max2_r" ) norm_tendon_length = 0.13374;
	else if( name=="glut_max3_r" ) norm_tendon_length = 0.1531;
	else if( name=="iliacus_r" ) norm_tendon_length = 0.10377;
	else if( name=="psoas_r" ) norm_tendon_length = 0.16455;
	else if( name=="quad_fem_r" ) norm_tendon_length = 0.025164;
	else if( name=="gem_r" ) norm_tendon_length = 0.041412;
	else if( name=="peri_r" ) norm_tendon_length = 0.11835;
	else if( name=="rect_fem_r" ) norm_tendon_length = 0.34369;
	else if( name=="vas_med_r" ) norm_tendon_length = 0.13779;
	else if( name=="vas_int_r" ) norm_tendon_length = 0.1499;
	else if( name=="vas_lat_r" ) norm_tendon_length = 0.17176;
	else if( name=="med_gas_r" ) norm_tendon_length = 0.38433;
	else if( name=="lat_gas_r" ) norm_tendon_length = 0.37406;
	else if( name=="soleus_r" ) norm_tendon_length = 0.24718;
	else if( name=="tib_post_r" ) norm_tendon_length = 0.30618;
	else if( name=="flex_dig_r" ) norm_tendon_length = 0.39849;
	else if( name=="flex_hal_r" ) norm_tendon_length = 0.37969;
	else if( name=="tib_ant_r" ) norm_tendon_length = 0.22124;
	else if( name=="per_brev_r" ) norm_tendon_length = 0.1585;
	else if( name=="per_long_r" ) norm_tendon_length = 0.34064;
	else if( name=="per_tert_r" ) norm_tendon_length = 0.099649;
	else if( name=="ext_dig_r" ) norm_tendon_length = 0.34414;
	else if( name=="ext_hal_r" ) norm_tendon_length = 0.30305;
	else if( name=="glut_med1_l" ) norm_tendon_length = 0.081673;
	else if( name=="glut_med2_l" ) norm_tendon_length = 0.055041;
	else if( name=="glut_med3_l" ) norm_tendon_length = 0.054844;
	else if( name=="glut_min1_l" ) norm_tendon_length = 0.017011;
	else if( name=="glut_min2_l" ) norm_tendon_length = 0.027745;
	else if( name=="glut_min3_l" ) norm_tendon_length = 0.054227;
	else if( name=="semimem_l" ) norm_tendon_length = 0.3958;
	else if( name=="semiten_l" ) norm_tendon_length = 0.27904;
	else if( name=="bifemlh_l" ) norm_tendon_length = 0.3589;
	else if( name=="bifemsh_l" ) norm_tendon_length = 0.096527;
	else if( name=="sar_l" ) norm_tendon_length = 0.1077;
	else if( name=="add_long_l" ) norm_tendon_length = 0.12192;
	else if( name=="add_brev_l" ) norm_tendon_length = 0.021815;
	else if( name=="add_mag1_l" ) norm_tendon_length = 0.066085;
	else if( name=="add_mag2_l" ) norm_tendon_length = 0.13573;
	else if( name=="add_mag3_l" ) norm_tendon_length = 0.28074;
	else if( name=="tfl_l" ) norm_tendon_length = 0.46862;
	else if( name=="pect_l" ) norm_tendon_length = 0.035865;
	else if( name=="grac_l" ) norm_tendon_length = 0.13845;
	else if( name=="glut_max1_l" ) norm_tendon_length = 0.12918;
	else if( name=="glut_max2_l" ) norm_tendon_length = 0.13324;
	else if( name=="glut_max3_l" ) norm_tendon_length = 0.15123;
	else if( name=="iliacus_l" ) norm_tendon_length = 0.10465;
	else if( name=="psoas_l" ) norm_tendon_length = 0.16577;
	else if( name=="quad_fem_l" ) norm_tendon_length = 0.025127;
	else if( name=="gem_l" ) norm_tendon_length = 0.04137;
	else if( name=="peri_l" ) norm_tendon_length = 0.11749;
	else if( name=="rect_fem_l" ) norm_tendon_length = 0.34832;
	else if( name=="vas_med_l" ) norm_tendon_length = 0.14179;
	else if( name=="vas_int_l" ) norm_tendon_length = 0.15312;
	else if( name=="vas_lat_l" ) norm_tendon_length = 0.17617;
	else if( name=="med_gas_l" ) norm_tendon_length = 0.38871;
	else if( name=="lat_gas_l" ) norm_tendon_length = 0.37655;
	else if( name=="soleus_l" ) norm_tendon_length = 0.25044;
	else if( name=="tib_post_l" ) norm_tendon_length = 0.30689;
	else if( name=="flex_dig_l" ) norm_tendon_length = 0.39916;
	else if( name=="flex_hal_l" ) norm_tendon_length = 0.38039;
	else if( name=="tib_ant_l" ) norm_tendon_length = 0.21818;
	else if( name=="per_brev_l" ) norm_tendon_length = 0.15862;
	else if( name=="per_long_l" ) norm_tendon_length = 0.34101;
	else if( name=="per_tert_l" ) norm_tendon_length = 0.0995;
	else if( name=="ext_dig_l" ) norm_tendon_length = 0.34089;
	else if( name=="ext_hal_l" ) norm_tendon_length = 0.30194;
	else if( name=="ercspn_r" ) norm_tendon_length = 0.032364;
	else if( name=="ercspn_l" ) norm_tendon_length = 0.032329;
	else if( name=="intobl_r" ) norm_tendon_length = 0.10689;
	else if( name=="intobl_l" ) norm_tendon_length = 0.10674;
	else if( name=="extobl_r" ) norm_tendon_length = 0.14704;
	else if( name=="extobl_l" ) norm_tendon_length = 0.14676;
	else if( name=="muscle_1" ) norm_tendon_length = 0.106712872;
	norm_tendon_length /= _optimalFiberLength;
	assert( norm_tendon_length >= 0 );
	if( _pennationAngle != 0.0 )
	{
		//double muscle_width = _optimalFiberLength * sin(_pennationAngle);
		//double J = ( norm_muscle_tendon_length - norm_tendon_length ) * _optimalFiberLength / muscle_width;
		double J = ( norm_muscle_tendon_length - norm_tendon_length ) / sin(_pennationAngle);
		ca = sqrt( J*J/(1+J*J) );
	}
	else // pennation is always zero for this muscle, so cos(alpha)=cos(_pennationAngle)=cos(0)=1
	{
		ca = 1.0;
	}
	normState[STATE_FIBER_LENGTH] = (norm_muscle_tendon_length - norm_tendon_length) / ca;
	_fiberLength = normState[STATE_FIBER_LENGTH] * _optimalFiberLength;
   _passiveForce = calcPassiveForce(normState[STATE_FIBER_LENGTH]);
	_activeForce = calcActiveForce(normState[STATE_FIBER_LENGTH]);
	normStateDeriv[STATE_FIBER_LENGTH] = 0.0; // this is the value sent to the integrator for derivative of fiber length
	double normFiberLengthDeriv = _speed / _optimalFiberLength * ca; // this is the actual fiber length derivative we use, only here
	double normForceVelocityForceValue = calcForceVelocity( normState[STATE_ACTIVATION], normFiberLengthDeriv );
   /* Un-normalize the muscle state derivatives and forces. */
   /* Note: Do not need to Un-Normalize activation dynamics equation since activation, deactivation parameters
	 specified in muscle file are now independent of time scale */
   _activationDeriv = normStateDeriv[STATE_ACTIVATION];
   _fiberLengthDeriv = normStateDeriv[STATE_FIBER_LENGTH] * Vmax;
	_passiveForce *= _maxIsometricForce;
	_activeForce *= normState[STATE_ACTIVATION] * _maxIsometricForce * normForceVelocityForceValue;
	_tendonForce = ( _activeForce + _passiveForce ) * ca;

	setForce(_tendonForce);
}

//_____________________________________________________________________________
/**
 * From cmg_dt.c - calc_tendon_force_dt
 *
 * CALC_TENDON_FORCE_DT: this routine calculates the force in tendon by finding
 * tendon strain and using it in an exponential function (JBME 2003 - Thelen)
 * FmaxTendonStrain - Function is parameterized by the tendon strain due to maximum isometric muscle force
 *     This should be specified as a dynamic parameter in the muscle file
 *
 * @param aNormTendonLength Normalized length of the tendon.
 * @return The force in the tendon.
 */
double Thelen2003MuscleMinusTendonForceLength::calcTendonForce(double aNormTendonLength) const
{
   double norm_resting_length = _tendonSlackLength / _optimalFiberLength;
   double tendon_strain =  (aNormTendonLength - norm_resting_length) / norm_resting_length;

	double KToe = 3;
	double ToeStrain = 0.609*_fmaxTendonStrain;
	double ToeForce = 0.333333;
	double klin = 1.712/_fmaxTendonStrain;

	double tendon_force;
	if (tendon_strain>ToeStrain)
		tendon_force = klin*(tendon_strain-ToeStrain)+ToeForce;
	else if (tendon_strain>0) 
		tendon_force = ToeForce*(exp(KToe*tendon_strain/ToeStrain)-1.0)/(exp(KToe)-1);
	else
		tendon_force=0.;

	// Add on a small stiffness so that tendon never truly goes slack for non-zero tendon lengths
	tendon_force+=0.001*(1.+tendon_strain);

   return tendon_force;
}

//_____________________________________________________________________________
/**
 * From gmc.dt.c - calc_passive_fiber_force_dt
 *
 * CALC_PASSIVE_FIBER_FORCE_DT: written by Darryl Thelen
 * this routine calculates the passive force in the muscle fibers using
 * an exponential-linear function instead of cubic splines.
 * It always returns a non-zero force for all muscle lengths
 * This equation is parameterized using the following dynamic parameters
 * which must be specified in the muscle file
 * Dynamic Parameters:
 *   FmaxMuscleStrain - passive muscle strain due to the application of 
 *                      maximum isometric muscle force
 *	 KshapePassive - exponential shape factor
 *
 *  The normalized force due to passive stretch is given by
 *  For L < (1+maxStrain)*Lo
 *		f/f0 = exp(ks*(L-1)/maxStrain) / exp(ks)
 *
 * @param aNormFiberLength Normalized length of the muscle fiber.
 * @return The passive force in the muscle fibers.
 */
double Thelen2003MuscleMinusTendonForceLength::calcPassiveForce(double aNormFiberLength) const
{
	double passive_force;

	if (aNormFiberLength>(1+_fmaxMuscleStrain)) { // Switch to a linear model at large forces
		double slope=(_kShapePassive/_fmaxMuscleStrain)*(exp(_kShapePassive*(1.0+_fmaxMuscleStrain-1.0)/_fmaxMuscleStrain)) / (exp(_kShapePassive));
		passive_force=1.0+slope*(aNormFiberLength-(1.0+_fmaxMuscleStrain));
	}
	else
		passive_force = (exp(_kShapePassive*(aNormFiberLength-1.0)/_fmaxMuscleStrain)) / (exp(_kShapePassive));

	return passive_force;
}

//_____________________________________________________________________________
/**
 * From gmc.dt.c - calc_active_force_dt
 *
 * CALC_ACTIVE_FORCE_DT: this routine calculates the active component of force
 * in the muscle fibers. It uses the current fiber length to interpolate the
 * active force-length curve - described by Gaussian curve as in Thelen, JBME 2003
 * *
 * @param aNormFiberLength Normalized length of the muscle fiber.
 * @return The active force in the muscle fibers.
 */
double Thelen2003MuscleMinusTendonForceLength::calcActiveForce(double aNormFiberLength) const
{
	double x=-(aNormFiberLength-1.)*(aNormFiberLength-1.)/_kShapeActive;
	return exp(x);
}

//_____________________________________________________________________________
/**
 * From gmc_dt.c - calc_norm_fiber_velocity_dt
 *
 * CALC_NORM_FIBER_VELOCITY_DT: written by Darryl Thelen
 * this routine calculates the normalized fiber velocity (scaled to Vmax) by inverting the
 * muscle force-velocity-activation relationship (Thelen, JBME 2003)
 * This equation is parameterized using the following dynamic parameters
 * which must be specified in the muscle file
 * Dynamic Parameters:
 *   damping - normalized passive damping in parallel with contractile element
 *   Af - velocity shape factor from Hill's equation
 *   Flen	- Maximum normalized force when muscle is lengthening
 *
 * @param aActivation Activation of the muscle.
 * @param aActiveForce Active force in the muscle fibers.
 * @param aVelocityDependentForce Force value that depends on fiber velocity.
 * @return The velocity of the muscle fibers.
 */
double Thelen2003MuscleMinusTendonForceLength::calcFiberVelocity(double aActivation, double aActiveForce, double aVelocityDependentForce) const
{
   double epsilon=1.e-6;

   // Don't allow zero activation
	if (aActivation<epsilon) 
		aActivation=epsilon;

	double Fa = aActivation*aActiveForce;
	double Fv = aVelocityDependentForce;

   double norm_fiber_velocity;
	if (Fv<Fa) {		// Muscle shortening
		if (Fv<0) {	// Extend the force-velocity curve for negative forces using linear extrapolation
			double F0=0;
			double b=Fa+F0/_af;
        	double fv0 = (F0-Fa)/(b+_damping);
			double F1=epsilon;
			b=Fa+F1/_af;
        	double fv1 = (F1-Fa)/(b+_damping);
			b = (F1-F0)/(fv1-fv0);
        	norm_fiber_velocity = fv0 + (Fv-F0)/b;
		}
		else {
			double b=Fa+Fv/_af;
			norm_fiber_velocity = (Fv-Fa)/(b+_damping);
		}
	}
	else if (Fv<(.95*Fa*_flen)) {
		double b=(2+2./_af)*(Fa*_flen-Fv)/(_flen-1.);
		norm_fiber_velocity = (Fv-Fa)/(b+_damping);
	}
	else {  // Extend the force-velocity curve for forces that exceed maximum using linear extrapolation
			double F0=.95*Fa*_flen;
			double b=(2+2./_af)*(Fa*_flen-F0)/(_flen-1.);
        	double fv0 = (F0-Fa)/(b+_damping);
			double F1=(.95+epsilon)*Fa*_flen;
			b=(2+2./_af)*(Fa*_flen-F1)/(_flen-1.);
        	double fv1 = (F1-Fa)/(b+_damping);
			b = (fv1-fv0)/(F1-F0);
        	norm_fiber_velocity = fv0 + b*(Fv-F0);
    }

    return norm_fiber_velocity;
}
//_____________________________________________________________________________
/**
 * Calculate normalized fiber force from the normalized fiber velocity using
 * the normalized force-velocity curve.
 *
 * @param aNormalizedVelocity Normalized velocity of the muscle fibers.
 * @return The value of the normalized force-velocity curve for the current normalized fiber velocity.
 */
double Thelen2003MuscleMinusTendonForceLength::calcForceVelocity(double aActivation, double aNormalizedVelocity) const
{
	double norm_force;

	if( aNormalizedVelocity < -1.0 ) {
		norm_force = ( aNormalizedVelocity + 1.0 ) / ( 1.0 + 1.0 / _af );
	}
	else if( aNormalizedVelocity < 0.0 ) {
		norm_force = ( aNormalizedVelocity + 1.0 ) / ( 1.0 - aNormalizedVelocity / _af );
	}
	else {
		double C = _flen - 1;
		double K = C / ( 2.0 + 2.0 / _af );
		double bound = K * ( 0.95*_flen - 1.0 ) / ( 1.0 - 0.95*_flen*_flen );
		if( aNormalizedVelocity < bound ) {
			norm_force = ( K + aNormalizedVelocity ) / ( K + _flen*aNormalizedVelocity );
		}
		else {
			double N = -0.95*0.95*_flen*_flen + 0.90*_flen;
			double D = 0.05*0.05*_flen*_flen / K;
			norm_force = ( aNormalizedVelocity*D - N ) / C;
		}
	}

	return norm_force;
}
/*double Thelen2003MuscleMinusTendonForceLength::calcForceVelocity(double aActivation, double aNormalizedVelocity) const
{
	double norm_force;

	double K = 0.25 + 0.75 * aActivation;
	if( aNormalizedVelocity <= 0.0 )
	{
		norm_force = ( aNormalizedVelocity + K ) / ( K - aNormalizedVelocity / _af );
	}
	else
	{
		double L = ( 2 + 2 / _af ) / ( _flen - 1 );
		norm_force = ( L * aNormalizedVelocity * _flen + K ) / ( L * aNormalizedVelocity + K );
	}

	if( norm_force < 0.0 ) norm_force = 0.0;

	return norm_force;
}*/
//_____________________________________________________________________________
/**
 * Get the stress in this actuator.  It is calculated as the force divided
 * by the maximum isometric force (which is proportional to its area).
 */
double Thelen2003MuscleMinusTendonForceLength::getStress() const
{
	return _force / _maxIsometricForce;
}

//_____________________________________________________________________________
/**
 * Find the force produced by an actuator (the musculotendon unit), assuming
 * static equilibrium. Using the total muscle-tendon length, it finds the
 * fiber and tendon lengths so that the forces in each match. This routine
 * takes pennation angle into account, so its definition of static equilibrium
 * is when tendon_force = fiber_force * cos(pennation_angle). This funcion
 * will modify the object's values for _length, _fiberLength, _activeForce, 
 * and _passiveForce.
 *
 * @param aActivation Activation of the muscle.
 * @return The isometric force in the muscle.
 */
double Thelen2003MuscleMinusTendonForceLength::
computeIsometricForce(double aActivation)
{
#define MAX_ITERATIONS 100
#define ERROR_LIMIT 0.01

   int i;
   double tendon_length, fiber_force, tmp_fiber_length, min_tendon_stiffness;
   double cos_factor, fiber_stiffness;
   double old_fiber_length, length_change, tendon_stiffness, percent;
   double error_force = 0.0, old_error_force, tendon_force, norm_tendon_length;
   
   // If the muscle has no fibers, then treat it as a ligament.
   if (_optimalFiberLength < ROUNDOFF_ERROR) {
		// ligaments should be a separate class, so _optimalFiberLength should
		// never be zero.
      return 0.0;
   }

	calcLengthAfterPathComputation();

   // Make first guess of fiber and tendon lengths. Make fiber length equal to
   // optimal_fiber_length so that you start in the middle of the active+passive
   // force-length curve. Muscle_width is the width, or thickness, of the
   // muscle-tendon unit. It is the shortest allowable fiber length because if
   // the muscle-tendon length is very short, the pennation angle will be 90
   // degrees and the fibers will be vertical (assuming the tendon is horizontal).
   // When this happens, the fibers are as long as the muscle is wide.
   // If the resting tendon length is zero, then set the fiber length equal to
   // the muscle tendon length / cosine_factor, and find its force directly.

   double muscle_width = _optimalFiberLength * sin(_pennationAngle);

   if (_tendonSlackLength < ROUNDOFF_ERROR) {
      tendon_length = 0.0;
      cos_factor = cos(atan(muscle_width / _length));
      _fiberLength = _length / cos_factor;

		_activeForce = calcActiveForce(_fiberLength / _optimalFiberLength) * aActivation;
		if (_activeForce < 0.0)
			_activeForce = 0.0;

		_passiveForce = calcPassiveForce(_fiberLength / _optimalFiberLength);
		if (_passiveForce < 0.0)
			_passiveForce = 0.0;

      return (_activeForce + _passiveForce) * _maxIsometricForce * cos_factor;
   } else if (_length < _tendonSlackLength) {
      _fiberLength = muscle_width;
      tendon_length = _length;
      return 0.0;
   } else {
      _fiberLength = _optimalFiberLength;
      cos_factor = cos(calcPennation(_fiberLength, _optimalFiberLength, _pennationAngle));  
      tendon_length = _length - _fiberLength * cos_factor;

      /* Check to make sure tendon is not shorter than its slack length. If it
       * is, set the length to its slack length and re-compute fiber length.
       */
      if (tendon_length < _tendonSlackLength) {
         tendon_length = _tendonSlackLength;
         cos_factor = cos(atan(muscle_width / (_length - tendon_length)));
         _fiberLength = (_length - tendon_length) / cos_factor;
         if (_fiberLength < muscle_width)
            _fiberLength = muscle_width;
      }
   }

   // Muscle-tendon force is found using an iterative method. First, you guess
   // the length of the muscle fibers and the length of the tendon, and
   // calculate their respective forces. If the forces match (are within
   // ERROR_LIMIT of each other), stop; else change the length guesses based
   // on the error and try again.
   for (i = 0; i < MAX_ITERATIONS; i++) {
		_activeForce = calcActiveForce(_fiberLength / _optimalFiberLength) * aActivation;
      if (_activeForce < 0.0)
         _activeForce = 0.0;

		_passiveForce = calcPassiveForce(_fiberLength / _optimalFiberLength);
      if (_passiveForce < 0.0)
         _passiveForce = 0.0;

      fiber_force = (_activeForce + _passiveForce) * _maxIsometricForce * cos_factor;

      norm_tendon_length = tendon_length / _optimalFiberLength;
      tendon_force = calcTendonForce(norm_tendon_length) * _maxIsometricForce;

      old_error_force = error_force;
 
      error_force = tendon_force - fiber_force;

      if (DABS(error_force) <= ERROR_LIMIT) // muscle-tendon force found!
         break;

      if (i == 0)
         old_error_force = error_force;

      if (DSIGN(error_force) != DSIGN(old_error_force)) {
         percent = DABS(error_force) / (DABS(error_force) + DABS(old_error_force));
         tmp_fiber_length = old_fiber_length;
         old_fiber_length = _fiberLength;
         _fiberLength += percent * (tmp_fiber_length - _fiberLength);
      } else {
         // Estimate the stiffnesses of the tendon and the fibers. If tendon
         // stiffness is too low, then the next length guess will overshoot
         // the equilibrium point. So we artificially raise it using the
         // normalized muscle force. (_activeForce+_passiveForce) is the
         // normalized force for the current fiber length, and we assume that
         // the equilibrium length is close to this current length. So we want
         // to get force = (_activeForce+_passiveForce) from the tendon as well.
         // We hope this will happen by setting the tendon stiffness to
         // (_activeForce+_passiveForce) times its maximum stiffness.
			double tendon_elastic_modulus = 1200.0;
			double tendon_max_stress = 32.0;

         tendon_stiffness = calcTendonForce(norm_tendon_length) *
				_maxIsometricForce / _tendonSlackLength;

         min_tendon_stiffness = (_activeForce + _passiveForce) *
	         tendon_elastic_modulus * _maxIsometricForce /
	         (tendon_max_stress * _tendonSlackLength);

         if (tendon_stiffness < min_tendon_stiffness)
            tendon_stiffness = min_tendon_stiffness;

         fiber_stiffness = _maxIsometricForce / _optimalFiberLength *
            (calcActiveForce(_fiberLength / _optimalFiberLength)  +
            calcPassiveForce(_fiberLength / _optimalFiberLength));

         // determine how much the fiber and tendon lengths have to
         // change to make the error_force zero. But don't let the
	      // length change exceed half the optimal fiber length because
	      // that's too big a change to make all at once.
         length_change = fabs(error_force/(fiber_stiffness / cos_factor + tendon_stiffness));

         if (fabs(length_change / _optimalFiberLength) > 0.5)
            length_change = 0.5 * _optimalFiberLength;

         // now change the fiber length depending on the sign of the error
         // and the sign of the fiber stiffness (which equals the sign of
         // the slope of the muscle's force-length curve).
         old_fiber_length = _fiberLength;

         if (error_force > 0.0)
            _fiberLength += length_change;
         else
            _fiberLength -= length_change;
      }

      cos_factor = cos(calcPennation(_fiberLength, _optimalFiberLength, _pennationAngle));
      tendon_length = _length - _fiberLength * cos_factor;

      // Check to make sure tendon is not shorter than its slack length. If it is,
      // set the length to its slack length and re-compute fiber length.
      if (tendon_length < _tendonSlackLength) {
         tendon_length = _tendonSlackLength;
         cos_factor = cos(atan(muscle_width / (_length - tendon_length)));
         _fiberLength = (_length - tendon_length) / cos_factor;
      }
   }

   return tendon_force;
}

//_____________________________________________________________________________
/**
 * Find the force produced by muscle under isokinetic conditions assuming
 * an infinitely stiff tendon.  That is, all the shortening velocity of the
 * actuator (the musculotendon unit) is assumed to be due to the shortening
 * of the muscle fibers alone.  This methods calls
 * computeIsometricForce() and so alters the internal member variables of this
 * muscle.
 *
 * Note that the current implementation approximates the effect of the
 * force-velocity curve.  It does not account for the shortening velocity
 * when it is solving for the equilibrium length of the muscle fibers.  And,
 * a generic representation of the force-velocity curve is used (as opposed
 * to the implicit force-velocity curve assumed by this model.
 *
 * @param aActivation Activation of the muscle.
 * @return Isokinetic force generated by the actuator.
 * @todo Reimplement this methods with more accurate representation of the
 * force-velocity curve.
 */
double Thelen2003MuscleMinusTendonForceLength::
computeIsokineticForceAssumingInfinitelyStiffTendon(double aActivation)
{
	double isometricForce = computeIsometricForce(aActivation);
	double normalizedLength = _fiberLength / _optimalFiberLength;
	double normalizedVelocity = - _speed / (_vmax * _optimalFiberLength);
	normalizedVelocity *= cos(_pennationAngle);
	double normalizedForceVelocity = evaluateForceLengthVelocityCurve(1.0,1.0,normalizedVelocity);
	
	return isometricForce * normalizedForceVelocity;
}
