NAME=$1
NELS=u1452@nelstor0.cbu.uib.no:/elixir-chr/nels/users/u1452/Projects/UiO_Dahl_Chromatin_2018/MadeleineFosslie_MF/200110_A00943.B.Project_Fosslie-Libs12-2020-01-06/
scp -i /knut/u1452@nelstor0.cbu.uib.no.key ${NELS}${NAME}.fastq.gz .
