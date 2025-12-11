function a_steering=steeringVec(N,theta)
a_steering=exp(1j*2*pi*1/2*(0:N-1)'*theta);