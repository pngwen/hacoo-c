{ pkgs }: {
	deps = [
   pkgs.cunit
		pkgs.clang_12
		pkgs.ccls
		pkgs.gdb
		pkgs.gnumake
	];
}