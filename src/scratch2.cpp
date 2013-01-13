
int rank[2], size[2], namelen, xranks[] = { 0 };

  char processor_name[MPI_MAX_PROCESSOR_NAME];

  MPI_Group group_world, group_garch;

  MPI_Comm comm_garch;

  int send_val, recv_val, send_val2, recv_val2;


  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank[0]);

  MPI_Comm_size(MPI_COMM_WORLD, &size[0]);

  MPI_Get_processor_name(processor_name, &namelen);


  MPI_Comm_group(MPI_COMM_WORLD, &group_world);

  MPI_Group_excl(group_world, 1, xranks, &group_garch);

  MPI_Comm_create(MPI_COMM_WORLD, group_garch, &comm_garch);

  printf("Hello world! I’m rank %d of %d on %s\n", rank[0], size[0], processor_name);

  if (rank[0]) {
	int rank1;
	MPI_Comm_rank(comm_garch, &rank1);
	cout << "rank1 = " << rank1 << endl;
  }
