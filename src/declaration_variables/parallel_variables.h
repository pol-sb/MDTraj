integer comm, taskid, numproc, ierror, message
integer blocksize, residu, first_particle, last_particle
integer num_interacts, inter_residu, inter_blocksize, first_inter, last_inter
integer, allocatable :: interact_list(:,:)
integer :: particle_range(2), interact_range(2)


