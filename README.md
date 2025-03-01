# PANSY Meteor Extension



https://github.com/user-attachments/assets/2e66e59c-f924-4e65-8808-d4a50f88fefd



The <a href="https://pansy.eps.s.u-tokyo.ac.jp/en/">PANSY</a> radar is a powerful MST radar operating in Antarctica at Syowa station. It is primarily used to study the troposphere and the lower thermosphere using radar echoes of turbulence. The radar is extremely capable in terms of sensitivity also for meteors. For the purpose of extending the scientific outcome of PANSY, a receiver capable of detecting meteor head echoes was shipped to Syowa station in late 2024, and became operation in January 2025. This will allow creating a continuous record of meteor head echoes during the remainder of the operational lifetime of the PANSY radar. This is a unique opportunity to observe the southern hemisphere meteor flux via the head echo technique, which is not likely to occur anything soon.

Software and engineering documantation related to the PANSY meteor head echo receiver. This is a software defined radio receiver that operates independently of the PANSY receiver, using the transmit pulse leakthrough as a phase and radar experiment sequence reference. 

<img width="681" alt="Screenshot 2025-02-06 at 22 45 44" src="https://github.com/user-attachments/assets/96143e2d-e476-4a23-aa50-909337f93215" />

<img width="913" alt="Screenshot 2025-02-04 at 22 42 37" src="https://github.com/user-attachments/assets/23685b9a-d38c-402e-ae0a-d1c42ffa0ac4" />

The video below shows the troposphere-stratosphere mode analyzed using the receiver measurements connected to the 1 kW 19 antenna mini PANSY in Shigaraki. The radar sees mainly airplanes and ground clutter.

https://github.com/user-attachments/assets/849bf707-b011-4c54-bbf1-1a141d7f4ec6

## Instructions

To run the python signal processing, use the standard environment, instead of the conda environment of SanDRA.  

> conda deactivate
