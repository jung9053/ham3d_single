    


    sprintf(fname,"./QuadData/periodic_check.plt",nn);
    fp = fopen(fname,"w");
    fprintf(fp,"VARIABLES=\"X\",\"Y\",\"Z\"\n");
    fprintf(fp,"ZONE ZONETYPE=FEBRICK N=%d E=%d DATAPACKING=BLOCK\n",g->nnodes,g->ncells); 
    fprintf(fp,"VARLOCATION = (1=NODAL, 2=NODAL, 3=NODAL)\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[3*i]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[3*i+1]);
  fprintf(fp,"\n");
  for(i=0;i<g->nnodes;i++)
    fprintf(fp,"%f\n",g->x[3*i+2]);
  for(i=0;i<g->ncells;i++)
    fprintf(fp,"%d %d %d %d %d %d %d %d\n",g->conn[8*i],g->conn[8*i+1],g->conn[8*i+2],g->conn[8*i+3],g->conn[8*i+4],g->conn[8*i+5],g->conn[8*i+6],g->conn[8*i+7]);

  fclose(fp);

