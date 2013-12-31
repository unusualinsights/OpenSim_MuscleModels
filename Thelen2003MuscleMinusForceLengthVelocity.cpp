// Thelen2003MuscleMinusForceLengthVelocity.cpp
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
#include "Thelen2003MuscleMinusForceLengthVelocity.h"
#include <OpenSim/Common/SimmMacros.h>
#include <OpenSim/Common/rdMath.h>
#include <OpenSim/Common/DebugUtilities.h>

//=============================================================================
// STATICS
//=============================================================================
using namespace std;
using namespace OpenSim;

const int Thelen2003MuscleMinusForceLengthVelocity::STATE_ACTIVATION = 0;
const int Thelen2003MuscleMinusForceLengthVelocity::STATE_FIBER_LENGTH = 1;

//=============================================================================
// CONSTRUCTOR(S) AND DESTRUCTOR
//=============================================================================
//_____________________________________________________________________________
/**
 * Default constructor.
 */
Thelen2003MuscleMinusForceLengthVelocity::Thelen2003MuscleMinusForceLengthVelocity() :
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
Thelen2003MuscleMinusForceLengthVelocity::~Thelen2003MuscleMinusForceLengthVelocity()
{
}

//_____________________________________________________________________________
/**
 * Copy constructor.
 *
 * @param aMuscle Thelen2003MuscleMinusForceLengthVelocity to be copied.
 */
Thelen2003MuscleMinusForceLengthVelocity::Thelen2003MuscleMinusForceLengthVelocity(const Thelen2003MuscleMinusForceLengthVelocity &aMuscle) :
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
 * @return Pointer to a copy of this Thelen2003MuscleMinusForceLengthVelocity.
 */
Object* Thelen2003MuscleMinusForceLengthVelocity::copy() const
{
	Thelen2003MuscleMinusForceLengthVelocity *musc = new Thelen2003MuscleMinusForceLengthVelocity(*this);
	return(musc);
}

//=============================================================================
// CONSTRUCTION METHODS
//=============================================================================
//_____________________________________________________________________________
/**
 * Copy data members from one Thelen2003MuscleMinusForceLengthVelocity to another.
 *
 * @param aMuscle Thelen2003MuscleMinusForceLengthVelocity to be copied.
 */
void Thelen2003MuscleMinusForceLengthVelocity::copyData(const Thelen2003MuscleMinusForceLengthVelocity &aMuscle)
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
	_tendonForceConstantOfIntegration = aMuscle._tendonForceConstantOfIntegration;
	_tendonForceConstantOfIntegrationCalculated = aMuscle._tendonForceConstantOfIntegrationCalculated;
	_initialActivation = aMuscle._initialActivation;
	_initialActivationSet = aMuscle._initialActivationSet;
}

//_____________________________________________________________________________
/**
 * Set the data members of this Thelen2003MuscleMinusForceLengthVelocity to their null values.
 */
void Thelen2003MuscleMinusForceLengthVelocity::setNull()
{
	setType("Thelen2003MuscleMinusForceLengthVelocity");

	setNumControls(1); setNumStates(2); setNumPseudoStates(0);
	bindControl(0, _excitation, "excitation");
	bindState(STATE_ACTIVATION, _activation, "activation");
	bindState(STATE_FIBER_LENGTH, _fiberLength, "fiber_length");

	_excitation = 0.0;
	_activation = 0.0;
	_fiberLength = 0.0;
	_tendonForceConstantOfIntegration = -1000.0;
	_tendonForceConstantOfIntegrationCalculated = false;
	_initialActivation = -1000.0;
	_initialActivationSet = false;
}

//_____________________________________________________________________________
/**
 * Connect properties to local pointers.
 */
void Thelen2003MuscleMinusForceLengthVelocity::setupProperties()
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
 * @param aModel model containing this Thelen2003MuscleMinusForceLengthVelocity.
 */
void Thelen2003MuscleMinusForceLengthVelocity::setup(Model* aModel)
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
 * a Thelen2003MuscleMinusForceLengthVelocity.
 *
 * @param aActuator Actuator to copy property values from.
 */
void Thelen2003MuscleMinusForceLengthVelocity::copyPropertyValues(AbstractActuator& aActuator)
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
Thelen2003MuscleMinusForceLengthVelocity& Thelen2003MuscleMinusForceLengthVelocity::operator=(const Thelen2003MuscleMinusForceLengthVelocity &aMuscle)
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
void Thelen2003MuscleMinusForceLengthVelocity::scale(const ScaleSet& aScaleSet)
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
void Thelen2003MuscleMinusForceLengthVelocity::postScale(const ScaleSet& aScaleSet)
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
double Thelen2003MuscleMinusForceLengthVelocity::getPennationAngle()
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
double Thelen2003MuscleMinusForceLengthVelocity::getFiberLength()
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
double Thelen2003MuscleMinusForceLengthVelocity::getNormalizedFiberLength()
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
double Thelen2003MuscleMinusForceLengthVelocity::getPassiveFiberForce()
{
	return _passiveForce;
}
//_____________________________________________________________________________
/**
 * Get the applied active force generated by the muscle fibers.
 * Thelen2003MuscleMinusForceLengthVelocity implements this function to apply the
 * active force that was actually computed instead of just subtracting the
 * passive fiber force from the total fiber force.
 *
 * @return Current active force of the muscle fibers.
 */
double Thelen2003MuscleMinusForceLengthVelocity::getAppliedActiveFiberForce()
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
void Thelen2003MuscleMinusForceLengthVelocity::computeStateDerivatives(double rDYDT[])
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
void Thelen2003MuscleMinusForceLengthVelocity::computeEquilibrium()
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
void Thelen2003MuscleMinusForceLengthVelocity::computeActuation()
{
	// Base Class (to calculate speed)
	AbstractMuscle::computeActuation();

   double normState[2], normStateDeriv[2], norm_tendon_length, ca;
   double norm_muscle_tendon_length, pennation_angle;

   /* Normalize the muscle states */
   normState[STATE_ACTIVATION] = _activation;
   normState[STATE_FIBER_LENGTH] = _fiberLength / _optimalFiberLength;

   if( !_initialActivationSet ) { // if this variable hasn't been set yet, set it
	   _initialActivation = _activation;
	   _initialActivationSet = true;
   }

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

	pennation_angle = AbstractMuscle::calcPennation(normState[STATE_FIBER_LENGTH], 1.0, _pennationAngle);
   ca = cos(pennation_angle);

   norm_muscle_tendon_length = getLength() / _optimalFiberLength;
   norm_tendon_length = norm_muscle_tendon_length - normState[STATE_FIBER_LENGTH] * ca;

   _tendonForce = calcTendonForce(norm_tendon_length);
   _passiveForce = calcPassiveForce(normState[STATE_FIBER_LENGTH]);
	_activeForce = calcActiveForce(normState[STATE_FIBER_LENGTH]);
	if( !_tendonForceConstantOfIntegrationCalculated ) { // if this variable hasn't been set yet, set it
		_tendonForceConstantOfIntegration = _tendonForce - ( _passiveForce + _activeForce * normState[STATE_ACTIVATION] ) * ca;
		_tendonForceConstantOfIntegrationCalculated = true;
	}
	_tendonForce -= _tendonForceConstantOfIntegration;

	// NOTE: SimmZajacMuscle has this check, but Darryl's muscle didn't seem to
	// if (_activeForce < 0.0) _activeForce = 0.0;
 
   /* If pennation equals 90 degrees, fiber length equals muscle width and fiber
    * velocity goes to zero.  Pennation will stay at 90 until tendon starts to
    * pull, then "stiff tendon" approximation is used to calculate approximate
    * fiber velocity.
    */
   if (EQUAL_WITHIN_ERROR(ca, 0.0))
   {
      if (EQUAL_WITHIN_ERROR(_tendonForce, 0.0))
      {
         normStateDeriv[STATE_FIBER_LENGTH] = 0.0;
			// ms->fiber_velocity = 0.0;
		}
		else 
		{
         double h = norm_muscle_tendon_length - _tendonSlackLength;
         double w = _optimalFiberLength * sin(_pennationAngle);
         double new_fiber_length = sqrt(h*h + w*w) / _optimalFiberLength;
			double new_pennation_angle = AbstractMuscle::calcPennation(new_fiber_length, 1.0, _pennationAngle);
         double new_ca = cos(new_pennation_angle);
         normStateDeriv[STATE_FIBER_LENGTH] = getSpeed() / (Vmax * new_ca);
		}
	}
   else
   {
      double velocity_dependent_force = _tendonForce / ca - _passiveForce;
      normStateDeriv[STATE_FIBER_LENGTH] = calcFiberVelocity(normState[STATE_ACTIVATION],normStateDeriv[STATE_ACTIVATION],normState[STATE_FIBER_LENGTH],_activeForce,pennation_angle);
   }

   /* Un-normalize the muscle state derivatives and forces. */
   /* Note: Do not need to Un-Normalize activation dynamics equation since activation, deactivation parameters
     specified in muscle file are now independent of time scale */
   _activationDeriv = normStateDeriv[STATE_ACTIVATION];
   _fiberLengthDeriv = normStateDeriv[STATE_FIBER_LENGTH] * Vmax;

	_tendonForce *= _maxIsometricForce;
	_passiveForce *= _maxIsometricForce;
	_activeForce *= normState[STATE_ACTIVATION] * _maxIsometricForce;
	//ms->tendon_length = norm_tendon_length*(*(ms->optimal_fiber_length));

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
double Thelen2003MuscleMinusForceLengthVelocity::calcTendonForce(double aNormTendonLength) const
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
 * Calculate the derivative of the tendon force with respect to the
 * *un-normalized* fiber length, even though the argument into this function
 * is the normalized fiber length.
 *
 * @param aNormTendonLength Normalized length of the tendon.
 * @return The derivative of the *un-normalized* force in the tendon.
 */
double Thelen2003MuscleMinusForceLengthVelocity::calcTendonForceDerivative(double aNormTendonLength) const
{
	double norm_resting_length = _tendonSlackLength / _optimalFiberLength;
	double tendon_strain =  (aNormTendonLength - norm_resting_length) / norm_resting_length;

	double KToe = 3;
	double ToeStrain = 0.609*_fmaxTendonStrain;
	double ToeForce = 0.333333;
	double klin = 1.712/_fmaxTendonStrain;

	double tendon_force_deriv;
	if (tendon_strain>ToeStrain)
		tendon_force_deriv = _maxIsometricForce * klin / _tendonSlackLength;
	else if (tendon_strain>0) {
		double coefficient = _maxIsometricForce*ToeForce*KToe / ( (exp(KToe)-1)*ToeStrain*_tendonSlackLength );
		double x = KToe / ToeStrain * tendon_strain;
		tendon_force_deriv = coefficient * exp(x);
	}
	else
		tendon_force_deriv = 0.;

	// Add on a small stiffness so that tendon never truly goes slack for non-zero tendon lengths
	tendon_force_deriv += _maxIsometricForce * 0.001 / _tendonSlackLength;

	return tendon_force_deriv;
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
double Thelen2003MuscleMinusForceLengthVelocity::calcPassiveForce(double aNormFiberLength) const
{
	double f = -1;
	string name = getName();
	if( name=="glut_med1_r" ) f = 1.6627;
	else if( name=="glut_med2_r" ) f = 3.7563;
	else if( name=="glut_med3_r" ) f = 5.1344;
	else if( name=="glut_min1_r" ) f = 4.0459;
	else if( name=="glut_min2_r" ) f = 4.4287;
	else if( name=="glut_min3_r" ) f = 2.4475;
	else if( name=="semimem_r" ) f = 0.53447;
	else if( name=="semiten_r" ) f = 3.1707;
	else if( name=="bifemlh_r" ) f = 1.1416;
	else if( name=="bifemsh_r" ) f = 8.5982;
	else if( name=="sar_r" ) f = 1.6796;
	else if( name=="add_long_r" ) f = 4.0417;
	else if( name=="add_brev_r" ) f = 3.8497;
	else if( name=="add_mag1_r" ) f = 0.74326;
	else if( name=="add_mag2_r" ) f = 0.36592;
	else if( name=="add_mag3_r" ) f = 0.60154;
	else if( name=="tfl_r" ) f = 23.2729;
	else if( name=="pect_r" ) f = 1.0801;
	else if( name=="grac_r" ) f = 1.2554;
	else if( name=="glut_max1_r" ) f = 0.24736;
	else if( name=="glut_max2_r" ) f = 0.36478;
	else if( name=="glut_max3_r" ) f = 0.31298;
	else if( name=="iliacus_r" ) f = 64.929;
	else if( name=="psoas_r" ) f = 42.1637;
	else if( name=="quad_fem_r" ) f = 5.4246;
	else if( name=="gem_r" ) f = 9.6482;
	else if( name=="peri_r" ) f = 2.8671;
	else if( name=="rect_fem_r" ) f = 50.9398;
	else if( name=="vas_med_r" ) f = 4.0656;
	else if( name=="vas_int_r" ) f = 4.5404;
	else if( name=="vas_lat_r" ) f = 6.0183;
	else if( name=="med_gas_r" ) f = 23.543;
	else if( name=="lat_gas_r" ) f = 13.2785;
	else if( name=="soleus_r" ) f = 83.2946;
	else if( name=="tib_post_r" ) f = 55.9983;
	else if( name=="flex_dig_r" ) f = 5.6784;
	else if( name=="flex_hal_r" ) f = 4.8762;
	else if( name=="tib_ant_r" ) f = 2.4899;
	else if( name=="per_brev_r" ) f = 8.0433;
	else if( name=="per_long_r" ) f = 17.0519;
	else if( name=="per_tert_r" ) f = 2.1423;
	else if( name=="ext_dig_r" ) f = 3.4674;
	else if( name=="ext_hal_r" ) f = 1.0886;
	else if( name=="glut_med1_l" ) f = 7.1441;
	else if( name=="glut_med2_l" ) f = 12.7896;
	else if( name=="glut_med3_l" ) f = 21.0929;
	else if( name=="glut_min1_l" ) f = 4.729;
	else if( name=="glut_min2_l" ) f = 6.7004;
	else if( name=="glut_min3_l" ) f = 6.6758;
	else if( name=="semimem_l" ) f = 5.6582;
	else if( name=="semiten_l" ) f = 12.7127;
	else if( name=="bifemlh_l" ) f = 11.3394;
	else if( name=="bifemsh_l" ) f = 9.6297;
	else if( name=="sar_l" ) f = 1.0538;
	else if( name=="add_long_l" ) f = 1.3484;
	else if( name=="add_brev_l" ) f = 2.9624;
	else if( name=="add_mag1_l" ) f = 0.94989;
	else if( name=="add_mag2_l" ) f = 0.75416;
	else if( name=="add_mag3_l" ) f = 2.8487;
	else if( name=="tfl_l" ) f = 3.1271;
	else if( name=="pect_l" ) f = 0.37046;
	else if( name=="grac_l" ) f = 1.3662;
	else if( name=="glut_max1_l" ) f = 1.0766;
	else if( name=="glut_max2_l" ) f = 2.0699;
	else if( name=="glut_max3_l" ) f = 3.4562;
	else if( name=="iliacus_l" ) f = 8.7153;
	else if( name=="psoas_l" ) f = 5.8562;
	else if( name=="quad_fem_l" ) f = 3.3916;
	else if( name=="gem_l" ) f = 4.0575;
	else if( name=="peri_l" ) f = 5.7916;
	else if( name=="rect_fem_l" ) f = 7.6156;
	else if( name=="vas_med_l" ) f = 3.0256;
	else if( name=="vas_int_l" ) f = 5.0735;
	else if( name=="vas_lat_l" ) f = 4.1731;
	else if( name=="med_gas_l" ) f = 14.2304;
	else if( name=="lat_gas_l" ) f = 6.5032;
	else if( name=="soleus_l" ) f = 39.6568;
	else if( name=="tib_post_l" ) f = 51.5364;
	else if( name=="flex_dig_l" ) f = 3.4854;
	else if( name=="flex_hal_l" ) f = 3.2185;
	else if( name=="tib_ant_l" ) f = 3.5948;
	else if( name=="per_brev_l" ) f = 6.9809;
	else if( name=="per_long_l" ) f = 17.8969;
	else if( name=="per_tert_l" ) f = 3.2202;
	else if( name=="ext_dig_l" ) f = 4.5123;
	else if( name=="ext_hal_l" ) f = 1.6449;
	else if( name=="ercspn_r" ) f = 60.3918;
	else if( name=="ercspn_l" ) f = 74.1622;
	else if( name=="intobl_r" ) f = 17.4517;
	else if( name=="intobl_l" ) f = 3.7896;
	else if( name=="extobl_r" ) f = 10.9826;
	else if( name=="extobl_l" ) f = 28.2757;
	//else if( name=="muscle_1" ) f = 2.781860259; // for block sine wave
	else if( name=="muscle_1" ) f = 223.6666; // for block passive right
	return f / _maxIsometricForce;
}

//_____________________________________________________________________________
/**
 * Calculates the derivative of the passive force in the muscle fibers
 * with respect to the *un-normalized* muscle fiber length, even though the
 * argument into this function is the normalized muscle fiber lenghth.
 *
 * @param aNormFiberLength Normalized length of the muscle fiber.
 * @return The derivative of the *un-normalized* passive force in the muscle fibers.
 */
double Thelen2003MuscleMinusForceLengthVelocity::calcPassiveForceDerivative(double aNormFiberLength) const
{
	return 0.0;
}

//_____________________________________________________________________________
/**
 * From gmc.dt.c - calc_active_force_dt
 *
 * CALC_ACTIVE_FORCE_DT: this routine calculates the active component of force
 * in the muscle fibers. It uses the current fiber length to interpolate the
 * active force-length curve - described by Gaussian curve as in Thelen, JBME 2003
 * 
 * @param aNormFiberLength Normalized length of the muscle fiber.
 * @return The active force in the muscle fibers.
 */
double Thelen2003MuscleMinusForceLengthVelocity::calcActiveForce(double aNormFiberLength) const
{
	double f = -1;
	string name = getName();
	if( name=="glut_med1_r" ) f = 353.223;
	else if( name=="glut_med2_r" ) f = 31.9618;
	else if( name=="glut_med3_r" ) f = 4.4409;
	else if( name=="glut_min1_r" ) f = 58.7498;
	else if( name=="glut_min2_r" ) f = 6.9491;
	else if( name=="glut_min3_r" ) f = 2.6226;
	else if( name=="semimem_r" ) f = 1.9605;
	else if( name=="semiten_r" ) f = 2.2128;
	else if( name=="bifemlh_r" ) f = 2.7137;
	else if( name=="bifemsh_r" ) f = 88.7345;
	else if( name=="sar_r" ) f = 1.9885;
	else if( name=="add_long_r" ) f = 256.6003;
	else if( name=="add_brev_r" ) f = 157.819;
	else if( name=="add_mag1_r" ) f = 10.044;
	else if( name=="add_mag2_r" ) f = 7.7204;
	else if( name=="add_mag3_r" ) f = 11.3145;
	else if( name=="tfl_r" ) f = 18.5203;
	else if( name=="pect_r" ) f = 31.5921;
	else if( name=="grac_r" ) f = 7.1633;
	else if( name=="glut_max1_r" ) f = 4.6044;
	else if( name=="glut_max2_r" ) f = 12.5684;
	else if( name=="glut_max3_r" ) f = 10.764;
	else if( name=="iliacus_r" ) f = 463.817;
	else if( name=="psoas_r" ) f = 410.2562;
	else if( name=="quad_fem_r" ) f = 20.8399;
	else if( name=="gem_r" ) f = 1.2121;
	else if( name=="peri_r" ) f = 0.01;//f = 0;
	else if( name=="rect_fem_r" ) f = 117.9426;
	else if( name=="vas_med_r" ) f = 39.3896;
	else if( name=="vas_int_r" ) f = 41.9482;
	else if( name=="vas_lat_r" ) f = 57.3132;
	else if( name=="med_gas_r" ) f = 7.0426;
	else if( name=="lat_gas_r" ) f = 0.01;//f = 0;
	else if( name=="soleus_r" ) f = 1037.5741;
	else if( name=="tib_post_r" ) f = 148.4525;
	else if( name=="flex_dig_r" ) f = 0.01;//f = 0;
	else if( name=="flex_hal_r" ) f = 0.01;//f = 0;
	else if( name=="tib_ant_r" ) f = 26.9358;
	else if( name=="per_brev_r" ) f = 4.4451;
	else if( name=="per_long_r" ) f = 54.7023;
	else if( name=="per_tert_r" ) f = 6.1769;
	else if( name=="ext_dig_r" ) f = 17.0308;
	else if( name=="ext_hal_r" ) f = 5.3847;
	else if( name=="glut_med1_l" ) f = 300.2423;
	else if( name=="glut_med2_l" ) f = 127.7018;
	else if( name=="glut_med3_l" ) f = 108.3856;
	else if( name=="glut_min1_l" ) f = 16.932;
	else if( name=="glut_min2_l" ) f = 16.216;
	else if( name=="glut_min3_l" ) f = 20.5284;
	else if( name=="semimem_l" ) f = 220.8211;
	else if( name=="semiten_l" ) f = 32.9821;
	else if( name=="bifemlh_l" ) f = 157.7847;
	else if( name=="bifemsh_l" ) f = 4.5815;
	else if( name=="sar_l" ) f = 5.579;
	else if( name=="add_long_l" ) f = 4.461;
	else if( name=="add_brev_l" ) f = 4.0869;
	else if( name=="add_mag1_l" ) f = 10.0905;
	else if( name=="add_mag2_l" ) f = 21.1676;
	else if( name=="add_mag3_l" ) f = 119.4812;
	else if( name=="tfl_l" ) f = 20.5171;
	else if( name=="pect_l" ) f = 1.9874;
	else if( name=="grac_l" ) f = 1.4218;
	else if( name=="glut_max1_l" ) f = 92.0897;
	else if( name=="glut_max2_l" ) f = 216.1572;
	else if( name=="glut_max3_l" ) f = 157.1061;
	else if( name=="iliacus_l" ) f = 108.075;
	else if( name=="psoas_l" ) f = 120.6198;
	else if( name=="quad_fem_l" ) f = 3.2062;
	else if( name=="gem_l" ) f = 2.3011;
	else if( name=="peri_l" ) f = 42.4631;
	else if( name=="rect_fem_l" ) f = 38.7168;
	else if( name=="vas_med_l" ) f = 125.6017;
	else if( name=="vas_int_l" ) f = 191.3993;
	else if( name=="vas_lat_l" ) f = 288.4524;
	else if( name=="med_gas_l" ) f = 0.01;//f = 0;
	else if( name=="lat_gas_l" ) f = 0.01;//f = 0;
	else if( name=="soleus_l" ) f = 20.4492;
	else if( name=="tib_post_l" ) f = 15.7415;
	else if( name=="flex_dig_l" ) f = 3.0792;
	else if( name=="flex_hal_l" ) f = 3.0686;
	else if( name=="tib_ant_l" ) f = 400.2874;
	else if( name=="per_brev_l" ) f = 6.7456;
	else if( name=="per_long_l" ) f = 12.8458;
	else if( name=="per_tert_l" ) f = 10.0841;
	else if( name=="ext_dig_l" ) f = 123.4483;
	else if( name=="ext_hal_l" ) f = 12.2444;
	else if( name=="ercspn_r" ) f = 501.0289;
	else if( name=="ercspn_l" ) f = 434.0711;
	else if( name=="intobl_r" ) f = 42.4286;
	else if( name=="intobl_l" ) f = 67.0842;
	else if( name=="extobl_r" ) f = 27.7182;
	else if( name=="extobl_l" ) f = 22.8879;
	//else if( name=="muscle_1" ) f = 0.001469318; // for block sine wave
	else if( name =="muscle_1" ) f = 0.1; // for block passive right
	if( _initialActivation == 0.0 ) return 0.0; // force should be 0 if activation 0
	return f / _maxIsometricForce / _initialActivation;
}

//_____________________________________________________________________________
/**
 * Calculates the derivative of the active force in the muscle fibers with
 * respect to the *normalized* muscle fiber length.
 * 
 * @param aNormFiberLength Normalized length of the muscle fiber.
 * @return The derivative of the *normalized* active force in the muscle fibers.
 */
double Thelen2003MuscleMinusForceLengthVelocity::calcActiveForceDerivative(double aNormFiberLength) const
{
	return 0.0;
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
 * @param aActiveForce Normalized active force in the muscle fibers.
 * @return The normalized velocity of the muscle fibers.
 */
double Thelen2003MuscleMinusForceLengthVelocity::calcFiberVelocity(double aActivation, double aActivationDeriv, double aNormFiberLength, double aActiveForce, double aPennationAngle) const
{
	double cosa = cos( aPennationAngle );
	double sina = sin( aPennationAngle );
	double cos2a = cosa * cosa;
	double sin2a = sina * sina;
	double cos3a = cos2a * cosa;
	double normMuscleTendonLength = _length / _optimalFiberLength;
	double normTendonLength = normMuscleTendonLength - aNormFiberLength * cosa;
	double tendonForceDeriv = calcTendonForceDerivative( normTendonLength );
	double term1 = tendonForceDeriv * getSpeed() / cosa;
	double term2 = -aActivationDeriv * _maxIsometricForce * aActiveForce;
	double term3 = _tendonForce * _maxIsometricForce * sin2a / ( aNormFiberLength * _optimalFiberLength * cos3a );
	double term4 = aActivation * _maxIsometricForce / _optimalFiberLength * calcActiveForceDerivative( aNormFiberLength );
	double term5 = tendonForceDeriv / cos2a;
	double term6 = calcPassiveForceDerivative( aNormFiberLength );
	double numerator = term1 + term2;
	double denominator = term3 + term4 + term5 + term6;
	double fiber_velocity = numerator / denominator;

	// Maximum contraction velocity is an activation-scaled value
	double Vmax = _vmax;
	if (aActivation<1.0)
		Vmax = _vmax0 + aActivation*(Vmax-_vmax0);
	Vmax = Vmax*_optimalFiberLength;
	double norm_fiber_velocity = fiber_velocity / Vmax;

    return norm_fiber_velocity;
}
//_____________________________________________________________________________
/**
 * Get the stress in this actuator.  It is calculated as the force divided
 * by the maximum isometric force (which is proportional to its area).
 */
double Thelen2003MuscleMinusForceLengthVelocity::getStress() const
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
double Thelen2003MuscleMinusForceLengthVelocity::
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
double Thelen2003MuscleMinusForceLengthVelocity::
computeIsokineticForceAssumingInfinitelyStiffTendon(double aActivation)
{
	double isometricForce = computeIsometricForce(aActivation);
	double normalizedLength = _fiberLength / _optimalFiberLength;
	double normalizedVelocity = - _speed / (_vmax * _optimalFiberLength);
	normalizedVelocity *= cos(_pennationAngle);
	double normalizedForceVelocity = evaluateForceLengthVelocityCurve(1.0,1.0,normalizedVelocity);
	
	return isometricForce * normalizedForceVelocity;
}
