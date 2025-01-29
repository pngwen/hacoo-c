{ pkgs }: {
	deps = [
   pkgs.llvmPackages_11.openmp
   pkgs.cunit
		pkgs.clang_12
		pkgs.ccls
		pkgs.gdb
		pkgs.gnumake
	];
}