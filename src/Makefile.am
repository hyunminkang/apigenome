## Additional flags to pass to aclocal when it is invoked automatically at
## make time. The ${ACLOCAL_FLAGS} variable is picked up from the environment
## to provide a way for the user to supply additional arguments.
ACLOCAL_AMFLAGS = -I m4 ${ACLOCAL_FLAGS}

#SUBDIRS = htslib

## Paths to include
#TABIX_VER = 0.2.5
#TABIX = tabix-$(TABIX_VER)
AM_CPPFLAGS = -pipe -D__STDC_LIMIT_MACROS -Wall -Wno-unused-local-typedefs -Wno-enum-compare -fpic -O2 -I../../htslib #-I../include/ #-I$(TABIX)/
#AM_CPPFLAGS = -fopenmp -I ${top_builddir}/../htslib/ -pipe -D__STDC_LIMIT_MACROS -Wall -Wno-unused-local-typedefs -Wno-enum-compare -fpic -O2 -L ${top_builddir}/../htslib #-I../include/ #-I$(TABIX)/
#AM_CPPFLAGS = -I ../htslib/ -I ../htslib/htslib -pipe -D__STDC_LIMIT_MACROS -Wall -Wno-unused-local-typedefs -Wno-enum-compare -fpic -g #-O2 #-I../include/ #-I$(TABIX)/

EIGEN_HEADER_FILES = Eigen/Cholesky Eigen/CholmodSupport Eigen/Core Eigen/Dense Eigen/Eigen Eigen/Eigenvalues Eigen/Geometry Eigen/Householder Eigen/IterativeLinearSolvers Eigen/Jacobi Eigen/LU Eigen/MetisSupport Eigen/OrderingMethods Eigen/PaStiXSupport Eigen/PardisoSupport Eigen/QR Eigen/QtAlignedMalloc Eigen/SPQRSupport Eigen/SVD Eigen/Sparse Eigen/SparseCholesky Eigen/SparseCore Eigen/SparseLU Eigen/SparseQR Eigen/StdDeque Eigen/StdList Eigen/StdVector Eigen/SuperLUSupport Eigen/UmfPackSupport Eigen/src/Cholesky/LDLT.h Eigen/src/Cholesky/LLT.h Eigen/src/Cholesky/LLT_LAPACKE.h Eigen/src/CholmodSupport/CholmodSupport.h Eigen/src/Core/Array.h Eigen/src/Core/ArrayBase.h Eigen/src/Core/ArrayWrapper.h Eigen/src/Core/Assign.h Eigen/src/Core/AssignEvaluator.h Eigen/src/Core/Assign_MKL.h Eigen/src/Core/BandMatrix.h Eigen/src/Core/Block.h Eigen/src/Core/BooleanRedux.h Eigen/src/Core/CommaInitializer.h Eigen/src/Core/ConditionEstimator.h Eigen/src/Core/CoreEvaluators.h Eigen/src/Core/CoreIterators.h Eigen/src/Core/CwiseBinaryOp.h Eigen/src/Core/CwiseNullaryOp.h Eigen/src/Core/CwiseTernaryOp.h Eigen/src/Core/CwiseUnaryOp.h Eigen/src/Core/CwiseUnaryView.h Eigen/src/Core/DenseBase.h Eigen/src/Core/DenseCoeffsBase.h Eigen/src/Core/DenseStorage.h Eigen/src/Core/Diagonal.h Eigen/src/Core/DiagonalMatrix.h Eigen/src/Core/DiagonalProduct.h Eigen/src/Core/Dot.h Eigen/src/Core/EigenBase.h Eigen/src/Core/ForceAlignedAccess.h Eigen/src/Core/Fuzzy.h Eigen/src/Core/GeneralProduct.h Eigen/src/Core/GenericPacketMath.h Eigen/src/Core/GlobalFunctions.h Eigen/src/Core/IO.h Eigen/src/Core/Inverse.h Eigen/src/Core/Map.h Eigen/src/Core/MapBase.h Eigen/src/Core/MathFunctions.h Eigen/src/Core/MathFunctionsImpl.h Eigen/src/Core/Matrix.h Eigen/src/Core/MatrixBase.h Eigen/src/Core/NestByValue.h Eigen/src/Core/NoAlias.h Eigen/src/Core/NumTraits.h Eigen/src/Core/PermutationMatrix.h Eigen/src/Core/PlainObjectBase.h Eigen/src/Core/Product.h Eigen/src/Core/ProductEvaluators.h Eigen/src/Core/Random.h Eigen/src/Core/Redux.h Eigen/src/Core/Ref.h Eigen/src/Core/Replicate.h Eigen/src/Core/ReturnByValue.h Eigen/src/Core/Reverse.h Eigen/src/Core/Select.h Eigen/src/Core/SelfAdjointView.h Eigen/src/Core/SelfCwiseBinaryOp.h Eigen/src/Core/Solve.h Eigen/src/Core/SolveTriangular.h Eigen/src/Core/SolverBase.h Eigen/src/Core/StableNorm.h Eigen/src/Core/Stride.h Eigen/src/Core/Swap.h Eigen/src/Core/Transpose.h Eigen/src/Core/Transpositions.h Eigen/src/Core/TriangularMatrix.h Eigen/src/Core/VectorBlock.h Eigen/src/Core/VectorwiseOp.h Eigen/src/Core/Visitor.h Eigen/src/Core/arch/AVX/Complex.h Eigen/src/Core/arch/AVX/MathFunctions.h Eigen/src/Core/arch/AVX/PacketMath.h Eigen/src/Core/arch/AVX/TypeCasting.h Eigen/src/Core/arch/AVX512/MathFunctions.h Eigen/src/Core/arch/AVX512/PacketMath.h Eigen/src/Core/arch/AltiVec/Complex.h Eigen/src/Core/arch/AltiVec/MathFunctions.h Eigen/src/Core/arch/AltiVec/PacketMath.h Eigen/src/Core/arch/CUDA/Complex.h Eigen/src/Core/arch/CUDA/Half.h Eigen/src/Core/arch/CUDA/MathFunctions.h Eigen/src/Core/arch/CUDA/PacketMath.h Eigen/src/Core/arch/CUDA/PacketMathHalf.h Eigen/src/Core/arch/CUDA/TypeCasting.h Eigen/src/Core/arch/Default/Settings.h Eigen/src/Core/arch/NEON/Complex.h Eigen/src/Core/arch/NEON/MathFunctions.h Eigen/src/Core/arch/NEON/PacketMath.h Eigen/src/Core/arch/SSE/Complex.h Eigen/src/Core/arch/SSE/MathFunctions.h Eigen/src/Core/arch/SSE/PacketMath.h Eigen/src/Core/arch/SSE/TypeCasting.h Eigen/src/Core/arch/ZVector/Complex.h Eigen/src/Core/arch/ZVector/MathFunctions.h Eigen/src/Core/arch/ZVector/PacketMath.h Eigen/src/Core/functors/AssignmentFunctors.h Eigen/src/Core/functors/BinaryFunctors.h Eigen/src/Core/functors/NullaryFunctors.h Eigen/src/Core/functors/StlFunctors.h Eigen/src/Core/functors/TernaryFunctors.h Eigen/src/Core/functors/UnaryFunctors.h Eigen/src/Core/products/GeneralBlockPanelKernel.h Eigen/src/Core/products/GeneralMatrixMatrix.h Eigen/src/Core/products/GeneralMatrixMatrixTriangular.h Eigen/src/Core/products/GeneralMatrixMatrixTriangular_BLAS.h Eigen/src/Core/products/GeneralMatrixMatrix_BLAS.h Eigen/src/Core/products/GeneralMatrixVector.h Eigen/src/Core/products/GeneralMatrixVector_BLAS.h Eigen/src/Core/products/Parallelizer.h Eigen/src/Core/products/SelfadjointMatrixMatrix.h Eigen/src/Core/products/SelfadjointMatrixMatrix_BLAS.h Eigen/src/Core/products/SelfadjointMatrixVector.h Eigen/src/Core/products/SelfadjointMatrixVector_BLAS.h Eigen/src/Core/products/SelfadjointProduct.h Eigen/src/Core/products/SelfadjointRank2Update.h Eigen/src/Core/products/TriangularMatrixMatrix.h Eigen/src/Core/products/TriangularMatrixMatrix_BLAS.h Eigen/src/Core/products/TriangularMatrixVector.h Eigen/src/Core/products/TriangularMatrixVector_BLAS.h Eigen/src/Core/products/TriangularSolverMatrix.h Eigen/src/Core/products/TriangularSolverMatrix_BLAS.h Eigen/src/Core/products/TriangularSolverVector.h Eigen/src/Core/util/BlasUtil.h Eigen/src/Core/util/Constants.h Eigen/src/Core/util/DisableStupidWarnings.h Eigen/src/Core/util/ForwardDeclarations.h Eigen/src/Core/util/MKL_support.h Eigen/src/Core/util/Macros.h Eigen/src/Core/util/Memory.h Eigen/src/Core/util/Meta.h Eigen/src/Core/util/NonMPL2.h Eigen/src/Core/util/ReenableStupidWarnings.h Eigen/src/Core/util/StaticAssert.h Eigen/src/Core/util/XprHelper.h Eigen/src/Eigenvalues/ComplexEigenSolver.h Eigen/src/Eigenvalues/ComplexSchur.h Eigen/src/Eigenvalues/ComplexSchur_LAPACKE.h Eigen/src/Eigenvalues/EigenSolver.h Eigen/src/Eigenvalues/GeneralizedEigenSolver.h Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h Eigen/src/Eigenvalues/HessenbergDecomposition.h Eigen/src/Eigenvalues/MatrixBaseEigenvalues.h Eigen/src/Eigenvalues/RealQZ.h Eigen/src/Eigenvalues/RealSchur.h Eigen/src/Eigenvalues/RealSchur_LAPACKE.h Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h Eigen/src/Eigenvalues/SelfAdjointEigenSolver_LAPACKE.h Eigen/src/Eigenvalues/Tridiagonalization.h Eigen/src/Geometry/AlignedBox.h Eigen/src/Geometry/AngleAxis.h Eigen/src/Geometry/EulerAngles.h Eigen/src/Geometry/Homogeneous.h Eigen/src/Geometry/Hyperplane.h Eigen/src/Geometry/OrthoMethods.h Eigen/src/Geometry/ParametrizedLine.h Eigen/src/Geometry/Quaternion.h Eigen/src/Geometry/Rotation2D.h Eigen/src/Geometry/RotationBase.h Eigen/src/Geometry/Scaling.h Eigen/src/Geometry/Transform.h Eigen/src/Geometry/Translation.h Eigen/src/Geometry/Umeyama.h Eigen/src/Geometry/arch/Geometry_SSE.h Eigen/src/Householder/BlockHouseholder.h Eigen/src/Householder/Householder.h Eigen/src/Householder/HouseholderSequence.h Eigen/src/IterativeLinearSolvers/BasicPreconditioners.h Eigen/src/IterativeLinearSolvers/BiCGSTAB.h Eigen/src/IterativeLinearSolvers/ConjugateGradient.h Eigen/src/IterativeLinearSolvers/IncompleteCholesky.h Eigen/src/IterativeLinearSolvers/IncompleteLUT.h Eigen/src/IterativeLinearSolvers/IterativeSolverBase.h Eigen/src/IterativeLinearSolvers/LeastSquareConjugateGradient.h Eigen/src/IterativeLinearSolvers/SolveWithGuess.h Eigen/src/Jacobi/Jacobi.h Eigen/src/LU/Determinant.h Eigen/src/LU/FullPivLU.h Eigen/src/LU/InverseImpl.h Eigen/src/LU/PartialPivLU.h Eigen/src/LU/PartialPivLU_LAPACKE.h Eigen/src/LU/arch/Inverse_SSE.h Eigen/src/MetisSupport/MetisSupport.h Eigen/src/OrderingMethods/Amd.h Eigen/src/OrderingMethods/Eigen_Colamd.h Eigen/src/OrderingMethods/Ordering.h Eigen/src/PaStiXSupport/PaStiXSupport.h Eigen/src/PardisoSupport/PardisoSupport.h Eigen/src/QR/ColPivHouseholderQR.h Eigen/src/QR/ColPivHouseholderQR_LAPACKE.h Eigen/src/QR/CompleteOrthogonalDecomposition.h Eigen/src/QR/FullPivHouseholderQR.h Eigen/src/QR/HouseholderQR.h Eigen/src/QR/HouseholderQR_LAPACKE.h Eigen/src/SPQRSupport/SuiteSparseQRSupport.h Eigen/src/SVD/BDCSVD.h Eigen/src/SVD/JacobiSVD.h Eigen/src/SVD/JacobiSVD_LAPACKE.h Eigen/src/SVD/SVDBase.h Eigen/src/SVD/UpperBidiagonalization.h Eigen/src/SparseCholesky/SimplicialCholesky.h Eigen/src/SparseCholesky/SimplicialCholesky_impl.h Eigen/src/SparseCore/AmbiVector.h Eigen/src/SparseCore/CompressedStorage.h Eigen/src/SparseCore/ConservativeSparseSparseProduct.h Eigen/src/SparseCore/MappedSparseMatrix.h Eigen/src/SparseCore/SparseAssign.h Eigen/src/SparseCore/SparseBlock.h Eigen/src/SparseCore/SparseColEtree.h Eigen/src/SparseCore/SparseCompressedBase.h Eigen/src/SparseCore/SparseCwiseBinaryOp.h Eigen/src/SparseCore/SparseCwiseUnaryOp.h Eigen/src/SparseCore/SparseDenseProduct.h Eigen/src/SparseCore/SparseDiagonalProduct.h Eigen/src/SparseCore/SparseDot.h Eigen/src/SparseCore/SparseFuzzy.h Eigen/src/SparseCore/SparseMap.h Eigen/src/SparseCore/SparseMatrix.h Eigen/src/SparseCore/SparseMatrixBase.h Eigen/src/SparseCore/SparsePermutation.h Eigen/src/SparseCore/SparseProduct.h Eigen/src/SparseCore/SparseRedux.h Eigen/src/SparseCore/SparseRef.h Eigen/src/SparseCore/SparseSelfAdjointView.h Eigen/src/SparseCore/SparseSolverBase.h Eigen/src/SparseCore/SparseSparseProductWithPruning.h Eigen/src/SparseCore/SparseTranspose.h Eigen/src/SparseCore/SparseTriangularView.h Eigen/src/SparseCore/SparseUtil.h Eigen/src/SparseCore/SparseVector.h Eigen/src/SparseCore/SparseView.h Eigen/src/SparseCore/TriangularSolver.h Eigen/src/SparseLU/SparseLU.h Eigen/src/SparseLU/SparseLUImpl.h Eigen/src/SparseLU/SparseLU_Memory.h Eigen/src/SparseLU/SparseLU_Structs.h Eigen/src/SparseLU/SparseLU_SupernodalMatrix.h Eigen/src/SparseLU/SparseLU_Utils.h Eigen/src/SparseLU/SparseLU_column_bmod.h Eigen/src/SparseLU/SparseLU_column_dfs.h Eigen/src/SparseLU/SparseLU_copy_to_ucol.h Eigen/src/SparseLU/SparseLU_gemm_kernel.h Eigen/src/SparseLU/SparseLU_heap_relax_snode.h Eigen/src/SparseLU/SparseLU_kernel_bmod.h Eigen/src/SparseLU/SparseLU_panel_bmod.h Eigen/src/SparseLU/SparseLU_panel_dfs.h Eigen/src/SparseLU/SparseLU_pivotL.h Eigen/src/SparseLU/SparseLU_pruneL.h Eigen/src/SparseLU/SparseLU_relax_snode.h Eigen/src/SparseQR/SparseQR.h Eigen/src/StlSupport/StdDeque.h Eigen/src/StlSupport/StdList.h Eigen/src/StlSupport/StdVector.h Eigen/src/StlSupport/details.h Eigen/src/SuperLUSupport/SuperLUSupport.h Eigen/src/UmfPackSupport/UmfPackSupport.h Eigen/src/misc/Image.h Eigen/src/misc/Kernel.h Eigen/src/misc/RealSvd2x2.h Eigen/src/misc/blas.h Eigen/src/misc/lapack.h Eigen/src/misc/lapacke.h Eigen/src/misc/lapacke_mangling.h Eigen/src/plugins/ArrayCwiseBinaryOps.h Eigen/src/plugins/ArrayCwiseUnaryOps.h Eigen/src/plugins/BlockMethods.h Eigen/src/plugins/CommonCwiseBinaryOps.h Eigen/src/plugins/CommonCwiseUnaryOps.h Eigen/src/plugins/MatrixCwiseBinaryOps.h Eigen/src/plugins/MatrixCwiseUnaryOps.h

noinst_HEADERS = hDist.h \
	bam_endian.h \
	boolParser.h \
	cdsStat.h ncStat.h \
	TypeConversion.h \
	Constant.h \
	kseq.h pException.h tabix.h \
	ksort.h pFile.h \
	genomeLoci.h genomePosition.h IO.h kstring.h \
	codonHelper.h fVcf.h \
	genomeScore.h khash.h \
	boolExprParser.h \
	knetfile.h \
	wFile.h \
	Constant.h MathConstant.h \
	bcf_variant_key.h bcf_filter_arg.h discrete_log_helper.h \
	var_dict.h \
	$(EIGEN_HEADER_FILES)

#bin_PROGRAMS = vcfast bgzip tabix tinyfast
#bin_PROGRAMS = vcfast tinyfast bcfast cramore
bin_PROGRAMS = cramore

#noinst_LIBRARIES = libtabix.a

## vcfast files
cramore_SOURCES = cramore.cpp \
	Error.cpp PhredHelper.cpp params.cpp commands.cpp \
	dropseq.cpp \
	hts_utils.cpp bcf_ordered_reader.cpp interval.cpp interval_tree.cpp utils.cpp genome_interval.cpp log_tool.cpp estimator.cpp \
	bcf_ordered_writer.cpp reference_sequence.cpp \
	bcf_chunked_reader.cpp genomeChunk.cpp \
	bcf_filtered_reader.cpp \
	bam_ordered_reader.cpp \
	nuclear_pedigree.cpp \
	genotype_concordance.cpp \
	filter.cpp \
	libsvm.cpp libsvm_gkm.cpp \
	cmd_bed_delta_svm_train.cpp \
	cmd_bed_matched_shuffle.cpp \
	cmd_bed_shuffle.cpp \
	cmd_vcf_mendel_dup_conc.cpp cmd_vcf_sample_summary.cpp cmd_vcf_squeeze.cpp cmd_vcf_delta_svm.cpp \
	cmd_vcf_extract.cpp \
	cmd_vcf_svd.cpp tsv_reader.cpp \
	cmd_cram_simul_contam.cpp cmd_cram_verify_pair_id.cpp cmd_cram_sparse_genotype.cpp cmd_cram_flagstat.cpp \
	cmd_cram_procap_detect.cpp cmd_cram_context_indel_analysis.cpp \
	cmd_sc_multinom_em.cpp cmd_sc_multinom_gibbs.cpp cmd_sc_map_stamps.cpp cmd_sc_kallisto_count.cpp \
	Sort.cpp MathVector.cpp MathGold.cpp MathGenMin.cpp MathMatrix.cpp \
	Random.cpp cmd_vcf_infer_ancestry.cpp ancestry_estimator.cpp \
	cmd_vcf_infer_isaf.cpp frequency_estimator.cpp cmd_vcf_update_sites.cpp \
	sam_filtered_reader.cpp sc_drop_seq.cpp cmd_cram_freemux.cpp \
	sam_ordered_writer.cpp \
	cmd_cram_simuxlet.cpp  \
	cmd_cram_mux_pileup.cpp \
	cmd_cram_dense_genotype.cpp joint_genotype_block_reader.cpp joint_genotype_block_record.cpp variant.cpp variant_manip.cpp augmented_bam_record.cpp sex_ploidy_map.cpp allele.cpp vntr.cpp \
	cmd_fasta_gc_content.cpp fastaGC.cpp \
	cmd_vcf_normalize_depth.cpp \
	cmd_bgen_to_vcf.cpp \
	cmd_cram_verify_bam.cpp contam_estimator.cpp \
	cmd_cram_demuxlet.cpp \
	cmd_cram_update_rg.cpp

#cramore_LDADD = ../../htslib/libhts.a -lpthread -llzma -lz -lbz2 -lcurl -lcrypto
#cramore_LDFLAGS = -all-static -Wl,-Bdynamic,-L../../htslib,-lhts,-Bstatic
cramore_LDADD = ../../htslib/libhts.a -lpthread -llzma -lz -lbz2 -lgomp -lcurl -lcrypto

## vcfast files
#bcfast_SOURCES = bcfast.cpp Error.cpp PhredHelper.cpp params.cpp bcf_ordered_reader.cpp interval.cpp interval_tree.cpp utils.cpp hts_utils.cpp genome_interval.cpp
#bcfast_LDADD = ../../htslib/libhts.a  -lpthread -llzma -lz -lbz2 -lgomp -lcurl -lcrypto

## vcfast files
#vcfast_SOURCES = vcfast.cpp Error.cpp PhredHelper.cpp params.cpp codonHelper.cpp genomePosition.cpp IO.cpp
#vcfast_LDADD = libtabix.a

#tinyfast_SOURCES = tinyfast.cpp Error.cpp params.cpp
#tinyfast_LDADD = libtabix.a

#bgzip_SOURCES = bgzip.c bedidx.c bgzf.c knetfile.c kstring.c

#tabix_SOURCES = tabix.c bedidx.c bgzf.c knetfile.c kstring.c index.c

#libtabix_a_SOURCES = bedidx.c bgzf.c knetfile.c kstring.c index.c
