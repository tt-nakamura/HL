OBJ = HL.o solvde.o ludcmp.o zbrent.o ion.o mixlen.o opacity.o eos.o fdjac.o

init: $(OBJ) HL_shoot.o init.o newt.o odeint.o
	g++ $(OBJ) HL_shoot.o init.o newt.o odeint.o
fig3: $(OBJ) fig3.o
	g++ $(OBJ) fig3.o
fig4: $(OBJ) fig4.o
	g++ $(OBJ) fig4.o
fig5: $(OBJ) fig5.o
	g++ $(OBJ) fig5.o
main_seq: main_seq.o
	g++ main_seq.o