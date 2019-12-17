function [flagG, flagP] = directed_partitioning(G, P, k)
	num_eigs = k+1;
	num_cluster = k;

  	[VG,lambdaG]=eigs(G,num_eigs,'sm');
  	[VP,lambdaP]=eigs(P,num_eigs,'sm');
 	VP=smoothVector(G, VP, 0.7, 10);
	VP = gram_schmidt(VP);
	[flagG,centerG]=kmeans(VG,num_cluster);
	[flagP, centerP]=kmeans(VP, num_cluster)
end
