NEURON {
SUFFIX Kv14_6States_temperature2
USEION k READ ek WRITE ik
RANGE g, gbar
}

UNITS { (mV) = (millivolt) }

PARAMETER {
Hkc=0.078666 (/degC)
Hki=0.1224 (/degC)
Hkn=-0.43947 (/degC)
Hri=0.56642 (/degC)
Hric_n=0.23032 (/degC)
Hrn=-0.14856 (/degC)
Hvc=0.501 (/degC)
Hvic_c=1.0 (/degC)
Hvic_n=-0.41335 (/degC)
Hvn_ic=1.0 (/degC)
Ri=8.6163 ()
Rn=0.1897 ()
Vc= -13.2731 (mV)
Zc= 1.5254 ()
kc= 0.05117 (/ms)
ki= 0.010217 (/ms)
kn= 0.019945 (/ms)
ric_c=1.0 ()
ric_n=137.5808 ()
vc_ic=10.0 ()
vic_c=22.8906 ()
vic_n=0.0036417 ()
vn_ic=0.15513 ()
gbar = 36     (millimho/cm2)
}

ASSIGNED {
v    (mV)
ek   (mV)
g    (millimho/cm2)
ik   (milliamp/cm2)
kco  (/ms)
koc  (/ms)
}

STATE { CS OS IC1 CIC1 IN INC1 }

LOCAL F, A,rc_ic,rn_ic,Hrn_ic

BREAKPOINT {
SOLVE kin METHOD sparse
g = gbar*OS
ik = g*(v - ek)*(1e-3)
}

INITIAL { 
F= 96.485(joule/mV)
A= 8.134(joule/degC)*(celsius+273.15)
rc_ic=ric_c
rn_ic=ric_n
Hrn_ic=Hric_n
SOLVE kin STEADYSTATE sparse 
}

KINETIC kin {
rates(v)
~ CS <-> OS     (kco, koc)
~ CIC1 <-> IC1  (vc_ic*kco,rc_ic*vc_ic*koc)
~ CS <-> CIC1   (ric_c*vic_c*ki*exp(Hki*(celsius-25))*exp(Hvic_c*(celsius-25)), vic_c*ki*Ri*exp(Hki*(celsius-25))*exp(Hri*(celsius-25))*exp(Hvic_c*(celsius-25)))
~ OS <-> IC1    (ki*exp(Hki*(celsius-25)), ki*Ri*exp(Hki*(celsius-25))*exp(Hri*(celsius-25)))
~ OS <-> IN     (kn*exp(Hkn*(celsius-25)), kn*Rn*exp(Hkn*(celsius-25))*exp(Hrn*(celsius-25)))
~ IC1 <-> INC1  (vn_ic*rn_ic*kn*exp(Hkn*(celsius-25))*exp(Hrn_ic*(celsius-25))*exp(Hvn_ic*(celsius-25)), vn_ic*kn*Rn*exp(Hkn*(celsius-25))*exp(Hrn*(celsius-25))*exp(Hvn_ic*(celsius-25)))
~ IN <-> INC1  (ric_n*vic_n*ki*exp(Hki*(celsius-25))*exp(Hric_n*(celsius-25))*exp(Hvic_n*(celsius-25)) , vic_n*ki*Ri*exp(Hki*(celsius-25))*exp(Hri*(celsius-25))*exp(Hvic_n*(celsius-25)))
CONSERVE CS+OS+IC1+CIC1+IN+INC1 = 1
}

PROCEDURE rates(v(millivolt)) {
kco=kc*exp(Zc*F/A*(v-Vc-Hvc*(celsius-25)))*exp(Hkc*(celsius-25))
koc=kc*exp(-Zc*F/A*(v-Vc-Hvc*(celsius-25)))*exp(Hkc*(celsius-25))
}




