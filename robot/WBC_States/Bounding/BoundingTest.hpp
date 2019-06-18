#ifndef BOUNDING_TEST_Cheetah
#define BOUNDING_TEST_Cheetah

#include <WBC_States/Test.hpp>
#include <Dynamics/Quadruped.h>

template <typename T> class StateProvider;

namespace BoundingPhase{
    constexpr int lift_up = 0;
    constexpr int posture_keeping = 1;
    constexpr int bounding = 2;
    constexpr int NUM_BOUNDING_PHASE = 3;
};


template <typename T>
class BoundingTest: public Test<T>{
    public:
        BoundingTest(FloatingBaseModel<T>* , const RobotType& );
        virtual ~BoundingTest();

        Vec3<T> _body_pos;
        Vec3<T> _body_vel;
        Vec3<T> _body_acc;

        Vec3<T> _body_ori_rpy;
        // This does not have any frame meaning. just derivation of rpy
        Vec3<T> _body_ang_vel; 


        std::string _folder_name;
    protected:
        virtual void _UpdateTestOneStep();
        virtual void _TestInitialization();
        virtual int _NextPhase(const int & phase);
        virtual void _UpdateExtraData(Cheetah_Extra_Data<T> * ext_data);

        void _SettingParameter();

        Controller<T>* _body_up_ctrl;
        Controller<T>* _posture_keeping;
        Controller<T>* _bounding;

        StateProvider<T>* _sp;
};

#endif
