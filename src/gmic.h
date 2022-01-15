/*
 #
 #  File        : gmic.h
 #                ( C++ header file )
 #
 #  Description : GREYC's Magic for Image Computing - G'MIC API file
 #                ( https://gmic.eu )
 #
 #  Note        : Include this file in your C++ source code, if you
 #                want to use the G'MIC interpreter in your own program.
 #
 #  Copyright   : David Tschumperlé
 #                ( https://tschumperle.users.greyc.fr/ )
 #
 #  Licenses    : This file is 'dual-licensed', you have to choose one
 #                of the two licenses below to apply.
 #
 #                CeCILL-C
 #                The CeCILL-C license is close to the GNU LGPL.
 #                ( http://cecill.info/licences/Licence_CeCILL-C_V1-en.html )
 #
 #            or  CeCILL v2.1
 #                The CeCILL license is compatible with the GNU GPL.
 #                ( http://cecill.info/licences/Licence_CeCILL_V2.1-en.html )
 #
 #  This software is governed either by the CeCILL or the CeCILL-C license
 #  under French law and abiding by the rules of distribution of free software.
 #  You can  use, modify and or redistribute the software under the terms of
 #  the CeCILL or CeCILL-C licenses as circulated by CEA, CNRS and INRIA
 #  at the following URL: "http://cecill.info".
 #
 #  As a counterpart to the access to the source code and  rights to copy,
 #  modify and redistribute granted by the license, users are provided only
 #  with a limited warranty  and the software's author,  the holder of the
 #  economic rights,  and the successive licensors  have only  limited
 #  liability.
 #
 #  In this respect, the user's attention is drawn to the risks associated
 #  with loading,  using,  modifying and/or developing or reproducing the
 #  software by the user in light of its specific status of free software,
 #  that may mean  that it is complicated to manipulate,  and  that  also
 #  therefore means  that it is reserved for developers  and  experienced
 #  professionals having in-depth computer knowledge. Users are therefore
 #  encouraged to load and test the software's suitability as regards their
 #  requirements in conditions enabling the security of their systems and/or
 #  data to be ensured and,  more generally, to use and operate it in the
 #  same conditions as regards security.
 #
 #  The fact that you are presently reading this means that you have had
 #  knowledge of the CeCILL and CeCILL-C licenses and that you accept its terms.
 #
*/

#ifndef gmic_version
#define gmic_version 301

#ifndef gmic_pixel_type
#define gmic_pixel_type float
#endif

#include <cstdio>
#include <cstring>

// Define gmic_uint64 type.
#ifndef gmic_uint64
#if cimg_OS==2
#define gmic_uint64 __int64
#else // #if cimg_OS==2
#if UINTPTR_MAX==0xffffffff || defined(__arm__) || defined(_M_ARM) || ((ULONG_MAX)==(UINT_MAX))
#define gmic_uint64 unsigned long long
#else
#define gmic_uint64 unsigned long
#endif // #if UINTPTR_MAX==0xffffffff || defined(__arm__) || defined(_M_ARM) || ((ULONG_MAX)==(UINT_MAX))
#endif // #if cimg_OS==2
#endif // #ifndef gmic_uint64

#ifndef gmic_build

// Define classes 'gmic_image<T>' and 'gmic_list<T>'.
//---------------------------------------------------
#ifndef cimg_version

#define gmic_image CImg
#define gmic_list CImgList

namespace cimg_library {

  // Class 'gmic_image<T>'.
  template<typename T> struct gmic_image {
    unsigned int _width;       // Number of image columns (dimension along the X-axis)
    unsigned int _height;      // Number of image lines (dimension along the Y-axis)
    unsigned int _depth;       // Number of image slices (dimension along the Z-axis)
    unsigned int _spectrum;    // Number of image channels (dimension along the C-axis)
    bool _is_shared;           // Tells if the data buffer has been allocated by another object
    T *_data;                  // Pointer to the first pixel value

    // Destructor.
    ~gmic_image() {
      if (!_is_shared) delete[] _data;
    }

    // Constructor.
    gmic_image():_width(0),_height(0),_depth(0),_spectrum(0),_is_shared(false),_data(0) { }

    // Allocate memory for specified image dimensions.
    gmic_image<T>& assign(const unsigned int size_x, const unsigned int size_y=1,
                          const unsigned int size_z=1, const unsigned int size_c=1);

    // Create image by copying existing buffer of t values.
    template<typename t>
    gmic_image<T>& assign(const t *const values, const unsigned int size_x, const unsigned int size_y=1,
                          const unsigned int size_z=1, const unsigned int size_c=1);

    gmic_image<T>& assign(const T *const values, const unsigned int size_x, const unsigned int size_y,
                          const unsigned int size_z, const unsigned int size_c, const bool is_shared);

    // Pixel access.
    operator T*() {
      return _data;
    }

    operator const T*() const {
      return _data;
    }

    T& operator()(const unsigned int x, const unsigned int y=0, const unsigned z=0, const unsigned c=0) {
      return _data[x + y*_width + z*_width*_height + c*_width*_height*_depth ];
    }

    const T& operator()(const unsigned int x, const unsigned int y=0, const unsigned z=0, const unsigned c=0) const {
      return _data[x + y*_width + z*_width*_height + c*_width*_height*_depth ];
    }

  };

  // Class 'gmic_list<T>'.
  template<typename T> struct gmic_list {
    unsigned int _width;           // Number of images in the list
    unsigned int _allocated_width; // Allocated items in the list (must be 2^N and >size)
    gmic_image<T> *_data;          // Pointer to the first image of the list

    // Destructor.
    ~gmic_list() {
      delete[] _data;
    }

    // Constructor.
    gmic_list():_width(0),_allocated_width(0),_data(0) {}

    // Allocate memory for specified list dimension.
    gmic_list<T>& assign(const unsigned int n);

    // Image access.
    operator gmic_image<T>*() {
      return _data;
    }

    operator const gmic_image<T>*() const {
      return _data;
    }

    gmic_image<T>& operator()(const unsigned int l) {
      return _data[l];
    }

    const gmic_image<T>& operator()(const unsigned int l) const {
      return _data[l];
    }

  };
}
#undef gmic_image
#undef gmic_list
#endif // #ifndef cimg_version

#else // #ifndef gmic_build

// Define private functions, used to compile libgmic.
//---------------------------------------------------

#ifndef cimg_verbosity
#define cimg_verbosity 1
#endif

#ifdef _MSC_VER
#pragma comment(linker,"/STACK:6291456")
#pragma inline_depth(2)
#endif

#include <locale>
#ifdef cimg_version
#error "[gmic] *** Error *** File 'CImg.h' has been already included (should have been done first in file 'gmic.h')."
#endif
#define cimg_plugin "gmic.cpp"
#define cimglist_plugin "gmic.cpp"

#ifdef cimg_use_abort
inline bool *gmic_abort_ptr(bool *const p_is_abort);
#define cimg_abort_init bool *const gmic_is_abort = ::gmic_abort_ptr(0); cimg::unused(gmic_is_abort)
#define cimg_abort_test if (*gmic_is_abort) throw CImgAbortException()
#endif

template<typename T>
inline double gmic_mp_dollar(const char *const str,
                             void *const p_list, const T& pixel_type);
#define cimg_mp_operator_dollar(str) \
  ::gmic_mp_dollar(str,&imglist,(T)0)

template<typename Ts, typename T>
inline double gmic_mp_get(Ts *const ptr, const unsigned int siz, const bool to_string, const char *const str,
                          void *const p_list, const T& pixel_type);
#define cimg_mp_func_get(ptr,siz,to_string,str) \
  return ::gmic_mp_get(ptr,siz,to_string,str,&mp.imglist,(T)0)

template<typename Ts, typename T>
inline double gmic_mp_set(Ts *const ptr, const unsigned int siz, const char *const str,
                          void *const p_list, const T& pixel_type);
#define cimg_mp_func_set(ptr,siz,str) \
  return ::gmic_mp_set(ptr,siz,str,&mp.imglist,(T)0)

template<typename T, typename Ts>
inline double gmic_mp_name(const unsigned int ind, Ts *const out_str, const unsigned int siz,
                           void *const p_list, const T& pixel_type);
#define cimg_mp_func_name(ind,out_str,siz) \
  return ::gmic_mp_name(ind,out_str,siz,&mp.imglist,(T)0)

template<typename T>
inline double gmic_mp_run(char *const str,
                          void *const p_list, const T& pixel_type);
#define cimg_mp_func_run(str) \
  return ::gmic_mp_run(str,&mp.imglist,(T)0)

template<typename Ts, typename T>
inline double gmic_mp_store(const Ts *const ptr, const unsigned int siz,
                            const unsigned int w, const unsigned int h, const unsigned int d, const unsigned int s,
                            const bool is_compressed, const char *const str,
                            void *const p_list, const T& pixel_type);
#define cimg_mp_func_store(ptr,siz,w,h,d,s,is_compressed,str) \
  return ::gmic_mp_store(ptr,siz,w,h,d,s,is_compressed,str,&mp.imglist,(T)0)

#ifndef cimg_display
#define cimg_display 0
#endif
#ifndef cimg_appname
#define cimg_appname "gmic"
#endif
#include "CImg.h"

#if cimg_OS==2
#include <process.h>
#include <psapi.h>

#elif cimg_OS==1
#include <cerrno>
#include <sys/resource.h>
#include <sys/syscall.h>
#include <signal.h>

#endif // #if cimg_OS==2

#endif // #ifndef gmic_build

// Define some special character codes used for replacement in double quoted strings.
const char gmic_dollar = 23, gmic_lbrace = 24, gmic_rbrace = 25, gmic_comma = 26, gmic_dquote = 28,
  gmic_store = 29; // <- this one is only used in variable names.

// Define main libgmic class 'gmic'.
//----------------------------------
#define gmic_image cimg_library::CImg
#define gmic_list cimg_library::CImgList

// Class 'gmic'.
struct gmic {

  // Destructor.
  ~gmic();

  // Constructors.
  gmic();

  template<typename T>
  gmic(const char *const commands_line,
       const char *const custom_commands=0,
       const bool include_stdlib=true,
       float *const p_progress=0, bool *const p_is_abort=0,
       const T& pixel_type=(T)0);

  template<typename T>
  gmic(const char *const commands_line,
       gmic_list<T>& images, gmic_list<char>& images_names, const char *const custom_commands=0,
       const bool include_stdlib=true, float *const p_progress=0, bool *const p_is_abort=0);

  // Run G'MIC pipeline on an already-constructed object.
  template<typename T>
  gmic& run(const char *const commands_line,
            float *const p_progress=0, bool *const p_is_abort=0,
            const T& pixel_type=(T)0);

  template<typename T>
  gmic& run(const char *const commands_line,
            gmic_list<T> &images, gmic_list<char> &images_names,
            float *const p_progress=0, bool *const p_is_abort=0);

  // These functions return (or init) G'MIC-specific paths.
  static const char* path_user(const char *const custom_path=0);
  static const char* path_rc(const char *const custom_path=0);
  static bool init_rc(const char *const custom_path=0);

  // Functions below should be considered as *private*, and should not be used in user's code.
  template<typename T>
  static bool search_sorted(const char *const str, const T& list, const unsigned int length, unsigned int &out_ind);
  template<typename T>
  static double mp_dollar(const char *const str,
                          void *const p_list, const T& pixel_type);
  template<typename Ts, typename T>
  static double mp_get(Ts *const ptr, const unsigned int siz, const bool to_string, const char *const str,
                       void *const p_list, const T& pixel_type);
  template<typename Ts, typename T>
  static double mp_set(Ts *const ptr, const unsigned int siz, const char *const str,
                       void *const p_list, const T& pixel_type);
  template<typename T, typename Ts>
  static double mp_name(const unsigned int ind, Ts *const out_str, const unsigned int siz,
                        void *const p_list, const T& pixel_type);
  template<typename T>
  static double mp_run(char *const str,
                       void *const p_list, const T& pixel_type);
  template<typename Ts, typename T>
  static double mp_store(const Ts *const ptr, const unsigned int siz,
                         const unsigned int w, const unsigned int h, const unsigned int d, const unsigned int s,
                         const bool is_compressed, const char *const str,
                         void *const p_list, const T& pixel_type);
  static bool get_debug_info(const char *const s, unsigned int &line_number, unsigned int &file_number);
  static int _levenshtein(const char *const s, const char *const t,
                          gmic_image<int>& d, const int i, const int j);
  static int levenshtein(const char *const s, const char *const t);
  static unsigned int hashcode(const char *const str, const bool is_variable);
  static bool command_has_arguments(const char *const command);
  static const char* basename(const char *const str);
  static char *strreplace_fw(char *const str);
  static char *strreplace_bw(char *const str);
  static unsigned int strescape(const char *const str, char *const res);
  static const gmic_image<char>& decompress_stdlib();
  static bool *abort_ptr(bool *const p_is_abort);

  template<typename T>
  void _gmic(const char *const commands_line,
             gmic_list<T>& images, gmic_list<char>& images_names,
             const char *const custom_commands, const bool include_stdlib,
             float *const p_progress, bool *const p_is_abort);

  gmic_image<char> get_variable(const char *const name,
                                const unsigned int *const variables_sizes=0,
                                const gmic_list<char> *const images_names=0) const;
  const char *set_variable(const char *const name, const char operation='=',
                           const char *const value=0, const double *const pvalue=0,
                           const unsigned int *const variables_sizes=0);
  const char *set_variable(const char *const name, const gmic_image<unsigned char>& value,
                           const unsigned int *const variables_sizes=0);

  gmic& add_commands(const char *const data_commands, const char *const commands_file=0,
                     const bool add_debug_info=false,
                     unsigned int *count_new=0, unsigned int *count_replaced=0,
                     bool *const is_entrypoint=0);
  gmic& add_commands(std::FILE *const file, const char *const filename=0,
                     const bool add_debug_info=false,
                     unsigned int *count_new=0, unsigned int *count_replaced=0,
                     bool *const is_entrypoint=0);

  gmic_image<char> callstack2string(const bool _is_debug=false) const;
  gmic_image<char> callstack2string(const gmic_image<unsigned int>& callstack_selection,
                                    const bool _is_debug=false) const;
  gmic_image<char> callstack2string(const gmic_image<unsigned int>* callstack_selection,
                                    const bool _is_debug=false) const;
  void pop_callstack(const unsigned int callstack_size);

  gmic_image<unsigned int> selection2cimg(const char *const string, const unsigned int indice_max,
                                          const gmic_list<char>& names, const char *const command,
                                          const bool is_selection=true);

  gmic_image<char>& selection2string(const gmic_image<unsigned int>& selection,
                                     const gmic_list<char>& images_names,
                                     const unsigned int display_selection,
                                     gmic_image<char>& res) const;

  gmic_list<char> commands_line_to_CImgList(const char *const commands_line);

  template<typename T>
  void _gmic_substitute_args(const char *const argument, const char *const argument0,
                             const char *const command, const char *const item,
                             const gmic_list<T>& images);

  gmic& print(const char *format, ...);
  gmic& error(const bool output_header, const char *format, ...);
  gmic& debug(const char *format, ...);

  template<typename T>
  gmic_image<char> substitute_item(const char *const source,
                                   gmic_list<T>& images, gmic_list<char>& images_names,
                                   gmic_list<T>& parent_images, gmic_list<char>& parent_images_names,
                                   const unsigned int *const variables_sizes,
                                   const gmic_image<unsigned int> *const command_selection,
                                   const bool is_image_expr);

  template<typename T>
  void wait_threads(void *const p_gmic_threads, const bool try_abort, const T& pixel_type);

  template<typename T>
  gmic& print(const gmic_list<T>& list, const gmic_image<unsigned int> *const callstack_selection,
              const char *format, ...);

  template<typename T>
  gmic& warn(const gmic_list<T>& list, const gmic_image<unsigned int> *const callstack_selection,
             const bool force_visible, const char *format, ...);

  template<typename T>
  gmic& error(const bool output_header, const gmic_list<T>& list,
              const gmic_image<unsigned int> *const callstack_selection,
              const char *const command, const char *format, ...);

  template<typename T>
  bool check_cond(const char *const expr, gmic_list<T>& images, const char *const command);

  template<typename T>
  gmic& debug(const gmic_list<T>& list, const char *format, ...);

  template<typename T>
  gmic& print_images(const gmic_list<T>& images,
                     const gmic_list<char>& images_names,
                     const gmic_image<unsigned int>& selection,
                     const bool is_header=true);
  template<typename T>
  gmic& display_images(const gmic_list<T>& images,
                       const gmic_list<char>& images_names,
                       const gmic_image<unsigned int>& selection,
                       unsigned int *const XYZ,
                       const bool exit_on_anykey);
  template<typename T>
  gmic& display_plots(const gmic_list<T>& images,
                      const gmic_list<char>& images_names,
                      const gmic_image<unsigned int>& selection,
                      const unsigned int plot_type, const unsigned int vertex_type,
                      const double xmin, const double xmax,
                      const double ymin, const double ymax,
                      const bool exit_on_anykey);
  template<typename T>
  gmic& display_objects3d(const gmic_list<T>& images,
                          const gmic_list<char>& images_names,
                          const gmic_image<unsigned int>& selection,
                          const gmic_image<unsigned char>& background3d,
                          const bool exit_on_anykey);
  template<typename T>
  gmic_image<T>& check_image(const gmic_list<T>& list, gmic_image<T>& img);
  template<typename T>
  const gmic_image<T>& check_image(const gmic_list<T>& list, const gmic_image<T>& img);

  template<typename T>
  gmic& remove_images(gmic_list<T>& images, gmic_list<char>& images_names,
                      const gmic_image<unsigned int>& selection,
                      const unsigned int start, const unsigned int end);

  template<typename T>
  gmic& _run(const gmic_list<char>& commands_line,
             gmic_list<T> &images, gmic_list<char> &images_names,
             float *const p_progress, bool *const p_is_abort);

  template<typename T>
  gmic& _run(const gmic_list<char>& commands_line, unsigned int& position,
             gmic_list<T>& images, gmic_list<char>&images_names,
             gmic_list<T>& parent_images, gmic_list<char>& parent_images_names,
             const unsigned int *const variables_sizes,
             bool *const is_noargs, const char *const parent_arguments,
             const gmic_image<unsigned int> *const command_selection);

  // Class variables.
  static const char *builtin_commands_names[];
  static gmic_image<int> builtin_commands_inds;
  static gmic_image<char> stdlib;
  static gmic_list<void*> list_p_is_abort;
  static bool is_display_available;

  gmic_list<char> *commands, *commands_names, *commands_has_arguments,
    *_variables, *_variables_names, **variables, **variables_names,
    commands_files, callstack;
  gmic_image<unsigned int> dowhiles, fordones, repeatdones;
  gmic_image<unsigned char> light3d;
  gmic_image<void*> display_windows;
  gmic_image<char> status;

  float focale3d, light3d_x, light3d_y, light3d_z, specular_lightness3d, specular_shininess3d, _progress, *progress;
  gmic_uint64 reference_time;
  unsigned int nb_dowhiles, nb_fordones, nb_repeatdones, nb_carriages_default, nb_carriages_stdout,
    debug_filename, debug_line, cimg_exception_mode;
  int verbosity,render3d, renderd3d, network_timeout;
  bool allow_entrypoint, is_change, is_debug, is_running, is_start, is_return, is_quit, is_double3d, is_debug_info,
    _is_abort, *is_abort, is_abort_thread;
  const char *starting_commands_line;
};

// Class 'gmic_exception'.
//------------------------
struct gmic_exception {
  gmic_image<char> _command, _message;

  // Constructors.
  gmic_exception() {}

  gmic_exception(const char *const command, const char *const message) {
    if (command) {
      _command.assign((unsigned int)std::strlen(command) + 1,1,1,1);
      std::strcpy(_command._data,command);
    }
    if (message) {
      _message.assign((unsigned int)std::strlen(message) + 1,1,1,1);
      std::strcpy(_message._data,message);
    }
  }

  // Return error message.
  const char *what() const { // Give the error message returned by the G'MIC interpreter.
    return _message._data?_message._data:"";
  }

  const char *command() const {
    return _command._data?_command._data:"";
  }
};

inline bool *gmic_abort_ptr(bool *const p_is_abort) {
  return gmic::abort_ptr(p_is_abort);
}

template<typename T>
inline double gmic_mp_dollar(const char *const str,
                             void *const p_list, const T& pixel_type) {
  return gmic::mp_dollar(str,p_list,pixel_type);
}

template<typename Ts, typename T>
inline double gmic_mp_get(Ts *const ptr, const unsigned int siz, const bool to_string, const char *const str,
                          void *const p_list, const T& pixel_type) {
  return gmic::mp_get(ptr,siz,to_string,str,p_list,pixel_type);
}

template<typename Ts, typename T>
inline double gmic_mp_set(Ts *const ptr, const unsigned int siz, const char *const str,
                          void *const p_list, const T& pixel_type) {
  return gmic::mp_set(ptr,siz,str,p_list,pixel_type);
}

template<typename T, typename Ts>
inline double gmic_mp_name(const unsigned int ind, Ts *const out_str, const unsigned int siz,
                           void *const p_list, const T& pixel_type) {
  return gmic::mp_name(ind,out_str,siz,p_list,pixel_type);
}

template<typename T>
inline double gmic_mp_run(char *const str,
                          void *const p_list, const T& pixel_type) {
  return gmic::mp_run(str,p_list,pixel_type);
}

template<typename Ts, typename T>
inline double gmic_mp_store(const Ts *const ptr, const unsigned int siz,
                            const unsigned int w, const unsigned int h, const unsigned int d, const unsigned int s,
                            const bool is_compressed, const char *const str,
                            void *const p_list, const T& pixel_type) {
  return gmic::mp_store(ptr,siz,w,h,d,s,is_compressed,str,p_list,pixel_type);
}

#endif // #ifndef gmic_version

// Local Variables:
// mode: c++
// End:
